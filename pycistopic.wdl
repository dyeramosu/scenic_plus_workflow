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

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:latest"
        Int preemptible = 0
        Int disk_space = 128
    }

    command <<<
    set -e 

    mkdir tmpdir
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
    models=run_cgs_models_mallet('/tmp/Mallet/bin/mallet',
                    cisTopic_obj,
                    n_topics=list(range(54, 63, 3)),
                    n_cpu=24,
                    n_iter=500, 
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    tmp_path='tmpdir/', #Use SCRATCH if many models or big data set
                    save_path=None)
    
    # save models
    pickle.dump(models, 
                open(os.path.join('pycistopic_output_wdl', 'PDAC_500_iter_LDA_54_63.pkl'), 'wb'))
    
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

    
    # use cell-topic probabilities to generate dimensionality reductions
    from pycisTopic.clust_vis import *
    run_umap(cisTopic_obj, target  = 'cell', scale=True)
   
    
    # topic binarization
    from pycisTopic.topic_binarization import *
    ## binarize topic-region distributions
    region_bin_topics = binarize_topics(cisTopic_obj, method='otsu', ntop=3000, plot=False)
    ## binarize cell-topic distributions
    binarized_cell_topic = binarize_topics(cisTopic_obj, target='cell', method='li', plot=False)

    ## save 
    with open(os.path.join('pycistopic_output_wdl' + 'binarized_cell_topic.pkl'), 'wb') as f:
        pickle.dump(binarized_cell_topic, f)
    with open(os.path.join('pycistopic_output_wdl' + 'binarized_topic_region.pkl'), 'wb') as f:
        pickle.dump(region_bin_topics, f)
    
    
    # identify DARs
    from pycisTopic.diff_features import *
    ## impute the region accessibility exploting the cell-topic and topic-region probabilities
    imputed_acc_obj = impute_accessibility(cisTopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
    ## log-normalize the imputed data
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
    ## identify highly variable regions (not mandatory but will speed up hypothesis testing for identifying DARs)
    variable_regions = find_highly_variable_features(normalized_imputed_acc_obj,
                                            min_disp = 0.05,
                                            min_mean = 0.0125, 
                                            max_mean = 3,
                                            max_disp = np.inf,
                                            n_bins=20, 
                                            n_top_features=None,
                                            plot=False)
    
    # Identify DARs between groups
    markers_dict= find_diff_features(cisTopic_obj, 
                                 imputed_acc_obj,
                                 variable='states',
                                 var_features=variable_regions,
                                 contrasts=None,
                                 adjpval_thr=0.05,
                                 log2fc_thr=np.log2(1.5),
                                 n_cpu=5,
                                 split_pattern = '-') 
    
    ## save
    with open(os.path.join('pycistopic_output_wdl' + 'Imputed_accessibility.pkl'), 'wb') as f:
        pickle.dump(imputed_acc_obj, f, protocol=4)
    with open(os.path.join('pycistopic_output_wdl' + 'DARs.pkl'), 'wb') as f:
        pickle.dump(markers_dict, f, protocol=4)
    
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