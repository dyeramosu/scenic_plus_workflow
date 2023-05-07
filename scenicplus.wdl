version 1.0

workflow SCENIC_PLUS {
    input {
        String output_dir # gbucket (make sure to add / to end)

        # runtime
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:1.0.0"
        Int preemptible = 0
        Int disk_space = 128
        
        # first task
        File adata_file
        File cistopic_file
        File menr_file # and second task

        # second task
        File toronto_tf_file
    }
    
    call create_scenic_plus_object {
        input:
            cpu = cpu,
            memory = memory,
            docker = docker,
            preemptible = preemptible,
            disk_space = disk_space,
            output_dir = output_dir,

            adata_file = adata_file,
            cistopic_file = cistopic_file,
            menr_file = menr_file
    }
    call run_scenic_plus {
        input:
            cpu = cpu,
            memory = memory,
            docker = docker,
            preemptible = preemptible,
            disk_space = disk_space,
            output_dir = output_dir,

            scplus_obj_task1 = create_scenic_plus_object.og_scenic_plus_object,
            toronto_tf_file = toronto_tf_file, 
            menr_file = menr_file
    }

    output { 
        File scenic_plus_output = run_scenic_plus.scenic_plus_object
    }
}

task create_scenic_plus_object {
    input {
        String output_dir

        File adata_file
        File cistopic_file
        File menr_file

        Int cpu 
        Int memory 
        String docker
        Int preemptible
        Int disk_space
    }

    command <<<
        set -e

        mkdir tmpdir
        mkdir scenic_plus_output_wdl

        python << CODE
        # imports
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import pyranges as pr
        import pickle
        import dill
        import os
        import sys
        import warnings

        # open files
        adata = sc.read_h5ad('~{adata_file}')
        cistopic_obj = dill.load(open('~{cistopic_file}', 'rb'))
        menr = dill.load(open('~{menr_file}', 'rb'))

        # create object 
        from scenicplus.scenicplus_class import create_SCENICPLUS_object
        scplus_obj = create_SCENICPLUS_object(
            GEX_anndata = adata.raw.to_adata(),
            cisTopic_obj = cistopic_obj,
            menr = menr,
            bc_transform_func = lambda x: f'{x}___cisTopic' # function to convert scATAC-seq barcodes to scRNA-seq ones
        )
        scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

        # save
        dill.dump(scplus_obj, open(os.path.join('scenic_plus_output_wdl', 'og_scplus_obj.pkl'), 'wb'), protocol=-1)

        CODE

        gsutil -m cp scenic_plus_output_wdl/og_scplus_obj.pkl ~{output_dir}
    >>>

    output {
        File og_scenic_plus_object = 'scenic_plus_output_wdl/og_scplus_obj.pkl'
    }

    runtime {
        docker: docker
        memory: memory + "G"
        bootDiskSizeGb: 12
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

task run_scenic_plus {
    input {
        String output_dir

        Int cpu 
        Int memory 
        String docker
        Int preemptible
        Int disk_space

        File scplus_obj_task1
        File toronto_tf_file # utoronto_human_tfs_v_1.01.txt
        File menr_file
    }

    command <<<
        set -e
        echo $(pwd)
        mkdir tmpdir
        mkdir scenic_plus_output_wdl

        wget -O tmpdir/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
        chmod +x tmpdir/bedToBigBed

        python << CODE
        # imports
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import pyranges as pr
        import pickle
        import dill
        import os
        import sys
        import warnings

        ########################### SCENIC PLUS OBJECT TO MUON FUNC ######################################
        """Convert from SCENIC+ class to MuData
        Given a SCENIC+ class this function will convert to a MuData object.
        """
        from scenicplus.scenicplus_class import SCENICPLUS
        from mudata import MuData, AnnData
        #from pycistarget._io import dict_motif_enrichment_results_to_mudata
        import numpy as np
        import pandas as pd
        from collections import OrderedDict
        from typing import Tuple

        def scenicplus_object_to_mudata(
            scplus_obj: SCENICPLUS,
            search_space_key: str = 'search_space',
            region_to_gene_key: str = 'region_to_gene',
            TF_to_gene_key: str = 'TF2G_adj',
            eRegulon_AUC_key: str = 'eRegulon_AUC',
            eRegulon_metadata_key: str = 'eRegulon_metadata') -> MuData:
            """
            Convert scplus_obj to MuData
            Parameters
            ----------
                scplus_obj: SCENICPLUS 
                    a scenicplus object
                search_space_key: str = 'search_space' 
                    key under which the search space is stored in .uns
                region_to_gene_key: str = 'region_to_gene' 
                    key under which the region to gene importances are stored in .uns
                TF_to_gene_key: str = 'TF2G_adj' 
                    key under which the TF to gene importances are stored in .uns
                eRegulon_AUC_key: str = 'eRegulon_AUC' 
                    key under which the eRegulon AUC values are stored in .uns
                eRegulon_metadata_key: str = 'eRegulon_metadata' 
                    key under which the eRegulon metadata is stored in .uns
            
            Returns
            -------
                Tuple[MuData, MuData]
                    Mudata with gene expression/region accessibility data and eRegulons and MuData containing motif enrichment results.
            
            """
            not_stored = set(scplus_obj.uns.keys()) - set([search_space_key, region_to_gene_key, TF_to_gene_key, eRegulon_AUC_key, eRegulon_metadata_key])
            print(
                f"Following items in scplus_obj.uns will not be stored, store them seperatly if you want to keep them.\n\t{', '.join(not_stored)}")
            mudata_constructor = {}
            
            #Construct ACC AnnData
            adata_ACC = AnnData(
                X = scplus_obj.X_ACC.T, dtype = np.int32,
                obs = pd.DataFrame(index = scplus_obj.cell_names),
                var = scplus_obj.metadata_regions.infer_objects())
            mudata_constructor['ACC'] = adata_ACC
            
            #Construct EXP AnnData
            adata_EXP = AnnData(
                X = scplus_obj.X_EXP, dtype = np.int32,
                obs = pd.DataFrame(index = scplus_obj.cell_names),
                var = scplus_obj.metadata_genes.infer_objects())
            mudata_constructor['EXP'] = adata_EXP

            #Construct eRegulon AUC AnnDatas
            adata_AUC_region = AnnData(
                X = np.array(scplus_obj.uns[eRegulon_AUC_key]['Region_based'], dtype = np.float32), dtype = np.float32,
                obs = pd.DataFrame(index = scplus_obj.uns[eRegulon_AUC_key]['Region_based'].index),
                var = pd.DataFrame(index = scplus_obj.uns[eRegulon_AUC_key]['Region_based'].columns))
            mudata_constructor['AUC_target_regions'] = adata_AUC_region
            adata_AUC_gene = AnnData(
                X = np.array(scplus_obj.uns[eRegulon_AUC_key]['Gene_based'], dtype = np.float32), dtype = np.float32,
                obs = pd.DataFrame(index = scplus_obj.uns[eRegulon_AUC_key]['Gene_based'].index),
                var = pd.DataFrame(index = scplus_obj.uns[eRegulon_AUC_key]['Gene_based'].columns))
            mudata_constructor['AUC_target_genes'] = adata_AUC_gene

            #construct uns
            uns = OrderedDict()
            uns['search_space'] = scplus_obj.uns[search_space_key].explode('Distance').infer_objects()
            uns['region_to_gene'] = scplus_obj.uns[region_to_gene_key].explode('Distance').infer_objects()
            uns['TF_to_gene'] = scplus_obj.uns[TF_to_gene_key].infer_objects()
            uns['eRegulon_metadata'] = scplus_obj.uns[eRegulon_metadata_key].infer_objects()

            mdata = MuData(
                mudata_constructor,
                obs = scplus_obj.metadata_cell.infer_objects(),
                obsm = {key: np.array(scplus_obj.dr_cell[key], dtype = np.float32) for key in scplus_obj.dr_cell.keys()},
                uns = uns)

            #mdata_menr = dict_motif_enrichment_results_to_mudata(scplus_obj.menr)

            return mdata#, mdata_menr

        ########################### SCENIC PLUS OBJECT TO MUON FUNC END ######################################

        # open scenic plus object from task1
        scplus_obj = dill.load(open('~{scplus_obj_task1}', 'rb'))
        # only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
        scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]

        # run scenic plus
        from scenicplus.wrappers.run_scenicplus import run_scenicplus
        try:
            run_scenicplus(
                scplus_obj = scplus_obj,
                variable = ['GEX_states'],
                species = 'hsapiens',
                assembly = 'hg38',
                tf_file = '~{toronto_tf_file}',
                save_path = 'scenic_plus_output_wdl/',
                biomart_host = 'http://sep2019.archive.ensembl.org/',
                upstream = [1000, 150000],
                downstream = [1000, 150000],
                calculate_TF_eGRN_correlation = False,
                calculate_DEGs_DARs = True,
                export_to_loom_file = True,
                export_to_UCSC_file = True,
                path_bedToBigBed = 'tmpdir/',
                n_cpu = 24,
                _temp_dir = os.path.abspath("tmpdir/"))
        except Exception as e:
            #in case of failure, still save the object
            #dill.dump(scplus_obj, open(os.path.join('scenic_plus_output_wdl', 'scplus_obj.pkl'), 'wb'), protocol=-1)
            #raise(e)
            x = scenicplus_object_to_mudata(scplus_obj, search_space_key='search_space', region_to_gene_key='region_to_gene', TF_to_gene_key='TF2G_adj', eRegulon_AUC_key='eRegulon_AUC', eRegulon_metadata_key='eRegulon_metadata')
        

        # save mudata object
        x.write(os.path.join('scenic_plus_output_wdl', 'mudata.h5mu'))

        # pickle mudata
        import mudata
        import pickle
        from scenicplus.scenicplus_class import SCENICPLUS

        mdata = mudata.read_h5mu(os.path.join('scenic_plus_output_wdl', 'mudata.h5mu'))
        menr = dill.load(open('~{menr_file}', 'rb'))

        scplus_obj2 = SCENICPLUS(
            X_ACC = mdata['ACC'].X.T,
            X_EXP = mdata['EXP'].X,
            metadata_regions = mdata['ACC'].var,
            metadata_genes = mdata['EXP'].var,
            metadata_cell = mdata.obs,
            menr = menr ,
            dr_cell = {k: dict(mdata.obsm)[k] for k in dict(mdata.obsm).keys() if k not in mdata.mod.keys()},
            dr_region = {},
            uns = mdata.uns)

        scplus_obj2.uns['eRegulon_AUC'] = {
            'Gene_based': mdata['AUC_target_genes'].to_df(),
            'Region_based': mdata['AUC_target_regions'].to_df()
        }

        # save pickled object
        dill.dump(scplus_obj2, open(os.path.join('scenic_plus_output_wdl', 'scplus_obj2.pkl'), 'wb'))

        CODE

        gsutil -m cp scenic_plus_output_wdl/scplus_obj2.pkl ~{output_dir}
        gsutil -m cp scenic_plus_output_wdl/mudata.h5mu ~{output_dir}
    >>>

    output {
        File scenic_plus_object = 'scenic_plus_output_wdl/scplus_obj2.pkl'
    }

    runtime {
        docker: docker
        memory: memory + "G"
        bootDiskSizeGb: 100
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}