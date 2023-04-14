version 1.0

workflow cisTopic {
    call create_pycistopic_object

    output {
        File pycistopic_output = create_pycistopic_object.create_pycistopic_output
    }
}

task create_pycistopic_object {
    input {
        String output_dir # gbucket
        File atac_data_og_file
        File adata_file

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:latest"
        Int preemptible = 0
        Int disk_space = 128
    }

    command <<<
    set -e 

    mkdir tmp
    mkdir pycistopic_output_wdl

    python << CODE
    # imports
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import pickle
    from pycisTopic.cistopic_class import *

    # read input h5ad files
    atac_data_og = sc.read_h5ad('~{atac_data_og_file}')
    adata_og = sc.read_h5ad('~{adata_file}')
    
    # subset atacseq and rnaseq data to cells that are in both
    overlap = list(set(atac_data_og.obs.index.to_list()) & set(adata_og.obs.index.to_list()))
    adata = adata_og[adata_og.obs.index.isin(overlap)]
    atac_data = atac_data_og[atac_data_og.obs.index.isin(overlap)]

    # create cisTopic object
    cell_data = adata.obs
    count_array = atac_data.X.toarray()
    count_df = pd.DataFrame(data=count_array, index=atac_data.obs_names, columns=atac_data.var_names)
    count_df = count_df.transpose()

    cisTopic_obj = create_cistopic_object(count_df)
    # add cell information 
    cisTopic_obj.add_cell_data(cell_data)

    # save cisTopic object
    pickle.dump(cisTopic_obj,
                open(os.path.join('pycistopic_output_wdl', 'cistopic_obj.pkl'), 'wb'))
    
    CODE

    tar -czvf pycistopic_output.tar.gz pycistopic_output_wdl
    gsutil rsync -r pycistopic_output_wdl ~{output_dir}

    >>>

    output {
        File create_pycistopic_output = 'pycistopic_output.tar.gz'
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