version 1.0

workflow SCENIC_PLUS {
    call scenic_plus

    output {
        File scenic_plus_output = scenic_plus.scenic_plus_output
    }
}

task scenic_plus {
    input {
        String output_dir # gbucket
        File tf_file
        String biomart_host

        File adata_file
        File cistopic_file
        File menr_file
        File bedToBigBed_file

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:latest"
        Int preemptible = 0
        Int disk_space = 128
    }

    command <<<
        set -e

        mkdir tmp
        mkdir PDAC_scenic_plus_output

        python << CODE
        #imports
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import pyranges as pr

        import pickle
        import dill
        import os
        import sys
        import warnings
        warnings.filterwarnings("ignore")

        _stderr = sys.stderr   
        null = open(os.devnull,'wb')

        #get input objects
        adata = sc.read_h5ad('~{adata_file}')
        cistopic_obj = dill.load(open('~{cistopic_file}', 'rb'))
        menr = dill.load(open('~{menr_file}', 'rb')) 

        #create scenic plus object
        from scenicplus.scenicplus_class import create_SCENICPLUS_object

        scplus_obj = create_SCENICPLUS_object(
            GEX_anndata = adata.raw.to_adata(),
            cisTopic_obj = cistopic_obj,
            menr = menr,
            bc_transform_func = lambda x: f'{x}___cisTopic'
        )
        scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
        scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]

        # run scenicplus
        from scenicplus.wrappers.run_scenicplus import run_scenicplus
        try:
            run_scenicplus(
                scplus_obj = scplus_obj,
                variable = ['GEX_states'],
                species = 'hsapiens',
                assembly = 'hg38',
                tf_file = '~{tf_file}',
                save_path = os.path.join('PDAC_scenic_plus_output'),
                biomart_host = '~{biomart_host}',
                upstream = [1000, 150000],
                downstream = [1000, 150000],
                calculate_TF_eGRN_correlation = True,
                calculate_DEGs_DARs = True,
                export_to_loom_file = True,
                export_to_UCSC_file = True,
                path_bedToBigBed = '~{bedToBigBed_file}',
                n_cpu = 24,
                _temp_dir = 'tmp')
        except Exception as e:
            #in case of failure, still save the object
            dill.dump(scplus_obj, open(os.path.join('PDAC_scenic_plus_output', 'scplus_obj.pkl'), 'wb'), protocol=-1)
            raise(e)

        CODE

        tar -czvf scenic_plus_output.tar.gz PDAC_scenic_plus_output
        gsutil rsync -r PDAC_scenic_plus_output ~{output_dir}
    >>>

    output {
        File scenic_plus_output = 'scenic_plus_output.tar.gz'
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