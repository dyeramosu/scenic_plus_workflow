version 1.0

workflow cisTopic {
    
    call run_pycistopic

    output {
        File pycistopic_output = run_pycistopic.pycistopic_object
    }
}

task run_pycistopic {
    input {
        String output_dir # gbucket
        File atac_data_og_file
        File adata_file
        File mallet_file

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
                open(os.path.join('pycistopic_output_wdl', 'cistopic_obj_pre_models.pkl'), 'wb'))
    
    # RUN MODELS
    models=run_cgs_models_mallet('~{mallet_file}',
                    cisTopic_obj,
                    n_topics=list(range(3, 75, 3)),
                    n_cpu=24,
                    n_iter=500, 
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    tmp_path='/tmp', #Use SCRATCH if many models or big data set
                    save_path=None)
    
    # save models
    pickle.dump(models, 
                open(os.path.join('pycistopic_output_wdl', 'PDAC_500_iter_LDA_3_75.pkl'), 'wb'))
    
    # evaluate models
    model=evaluate_models(models,
                     select_model=None, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False)
    
    # add model to cistopic object
    cisTopic_obj.add_LDA_model(model)

    # save cisTopic object
    pickle.dump(cisTopic_obj,
                open(os.path.join('pycistopic_output_wdl', 'cistopic_obj.pkl'), 'wb'))
   
    CODE

    tar -czvf pycistopic_output.tar.gz pycistopic_output_wdl
    gsutil rsync -r pycistopic_output_wdl ~{output_dir}

    >>>

    output {
        File pycistopic_object = 'pycistopic_output.tar.gz'
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