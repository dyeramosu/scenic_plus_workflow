version 1.0

workflow SCENIC_PLUS {
    input {
        String output_dir # gbucket (make sure to add / to end)

        # runtime
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:latest"
        Int preemptible = 0
        Int disk_space = 128
        
        # first task
        File adata_file
        File cistopic_file
        File menr_file

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
            toronto_tf_file = toronto_tf_file
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
    }

    command <<<
        set -e
        mkdir tmpdir
        echo pwd
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
                calculate_TF_eGRN_correlation = True,
                calculate_DEGs_DARs = True,
                export_to_loom_file = True,
                export_to_UCSC_file = True,
                path_bedToBigBed = 'tmpdir/',
                n_cpu = 24,
                _temp_dir = os.path.abspath("tmpdir/"))
        except Exception as e:
            #in case of failure, still save the object
            dill.dump(scplus_obj, open(os.path.join('scenic_plus_output_wdl', 'scplus_obj.pkl'), 'wb'), protocol=-1)
            raise(e)

        CODE

        gsutil -m cp scenic_plus_output_wdl/scplus_obj.pkl ~{output_dir}
    >>>

    output {
        File scenic_plus_object = 'scenic_plus_output_wdl/scplus_obj.pkl'
    }

    runtime {
        docker: docker
        memory: memory + "G"
        bootDiskSizeGb: 50
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}