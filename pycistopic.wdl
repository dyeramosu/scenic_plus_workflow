version 1.0

workflow cisTopic {

    input {
        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:1.0.0"
        Int preemptible = 0
        Int disk_space = 128

        String output_dir # gbucket (make sure to add / to end)

        File atac_data_og_file # barcoding-multiome-atac_for-Deepika.h5ad
        File adata_file # scenicplus_modified_gem.h5ad
    }
    
    call create_pycistopic_object {
        input:
            cpu = cpu,
            memory = memory,
            docker = docker,
            preemptible = preemptible, 
            disk_space = disk_space,
            output_dir = output_dir,

            atac_data_og_file = atac_data_og_file,
            adata_file = adata_file
    }
    call run_models_LDA {
        input:
            cpu = cpu,
            memory = memory,
            docker = docker,
            preemptible = preemptible, 
            disk_space = disk_space,
            output_dir = output_dir,

            cistopic_obj_task1 = create_pycistopic_object.og_cistopic_obj
    }
    call topic_binarization {
        input:
            cpu = cpu,
            memory = memory,
            docker = docker,
            preemptible = preemptible, 
            disk_space = disk_space,
            output_dir = output_dir,

            cistopic_obj_task2 = run_models_LDA.cistopic_obj_plus_model
    }

    output {
        File pycistopic_output = topic_binarization.DARs
    }
}

task create_pycistopic_object {
    input {
        File atac_data_og_file
        File adata_file
        
        Int cpu
        Int memory
        String docker
        Int preemptible
        Int disk_space
        String output_dir # gbucket
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
        atac_data = sc.read_h5ad('~{atac_data_og_file}')
        adata = sc.read_h5ad('~{adata_file}')
        
        # subset atacseq and rnaseq data to cells that are in both
        #overlap = list(set(atac_data_og.obs.index.to_list()) & set(adata_og.obs.index.to_list()))
        #adata = adata_og[adata_og.obs.index.isin(overlap)]
        #atac_data = atac_data_og[atac_data_og.obs.index.isin(overlap)]

        # create cisTopic object
        cell_data = adata.obs
        count_array = atac_data.X.toarray()
        count_df = pd.DataFrame(data=count_array, index=atac_data.obs_names, columns=atac_data.var_names)
        count_df = count_df.transpose()

        og_cistopic_obj = create_cistopic_object(count_df)
        # add cell information 
        og_cistopic_obj.add_cell_data(cell_data)

        # save cisTopic object
        pickle.dump(og_cistopic_obj,
                    open(os.path.join('pycistopic_output_wdl', 'og_cistopic_obj.pkl'), 'wb'))

        CODE

        gsutil -m cp pycistopic_output_wdl/og_cistopic_obj.pkl ~{output_dir} 
    >>>

    output {
        File og_cistopic_obj = 'pycistopic_output_wdl/og_cistopic_obj.pkl'
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

task run_models_LDA {
    input {
        File cistopic_obj_task1

        Int cpu
        Int memory
        String docker
        Int preemptible
        Int disk_space
        String output_dir
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

        # get object
        cistopic_obj = pickle.load(open('~{cistopic_obj_task1}', 'rb'))

        # RUN MODELS
        models=run_cgs_models_mallet('/tmp/Mallet/bin/mallet',
                        cistopic_obj,
                        n_topics=list(range(3, 60, 3)),
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
                    open(os.path.join('pycistopic_output_wdl', 'PDAC_500_iter_LDA_3_60.pkl'), 'wb'))
        
        # evaluate models
        model=evaluate_models(models,
                        select_model=None, 
                        return_model=True, 
                        metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                        plot_metrics=False)
        
        # add model to cistopic object
        cistopic_obj.add_LDA_model(model)

        # use cell-topic probabilities to generate dimensionality reductions
        from pycisTopic.clust_vis import *
        run_umap(cistopic_obj, target  = 'cell', scale=True)

        # save cisTopic object
        pickle.dump(cistopic_obj,
                    open(os.path.join('pycistopic_output_wdl', 'cistopic_obj_plus_model.pkl'), 'wb'))
        
        CODE

        gsutil -m cp pycistopic_output_wdl/PDAC_500_iter_LDA_3_60.pkl ~{output_dir}
        gsutil -m cp pycistopic_output_wdl/cistopic_obj_plus_model.pkl ~{output_dir}
    >>>

    output {
        File cistopic_obj_plus_model = 'pycistopic_output_wdl/cistopic_obj_plus_model.pkl'
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

task topic_binarization {
    input {
        File cistopic_obj_task2

        Int cpu
        Int memory
        String docker
        Int preemptible
        Int disk_space
        String output_dir
    }

    command  <<<
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

        # get object
        cistopic_obj = pickle.load(open('~{cistopic_obj_task2}', 'rb'))   
        
        # topic binarization
        from pycisTopic.topic_binarization import *
        ## binarize topic-region distributions
        region_bin_topics = binarize_topics(cistopic_obj, method='otsu', ntop=3000, plot=False)
        ## binarize cell-topic distributions
        binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=False)

        ## save 
        with open(os.path.join('pycistopic_output_wdl' , 'binarized_cell_topic.pkl'), 'wb') as f:
            pickle.dump(binarized_cell_topic, f)
        with open(os.path.join('pycistopic_output_wdl' , 'binarized_topic_region.pkl'), 'wb') as f:
            pickle.dump(region_bin_topics, f)
        
        # identify DARs
        from pycisTopic.diff_features import *
        ## impute the region accessibility exploting the cell-topic and topic-region probabilities
        imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
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
        markers_dict= find_diff_features(cistopic_obj, 
                                    imputed_acc_obj,
                                    variable='states',
                                    var_features=variable_regions,
                                    contrasts=None,
                                    adjpval_thr=0.05,
                                    log2fc_thr=np.log2(1.5),
                                    n_cpu=5,
                                    split_pattern = '-') 
        
        ## save
        with open(os.path.join('pycistopic_output_wdl' , 'Imputed_accessibility.pkl'), 'wb') as f:
            pickle.dump(imputed_acc_obj, f, protocol=4)
        with open(os.path.join('pycistopic_output_wdl' , 'DARs.pkl'), 'wb') as f:
            pickle.dump(markers_dict, f, protocol=4)
        
        CODE

        gsutil -m cp pycistopic_output_wdl/Imputed_accessibility.pkl ~{output_dir}
        gsutil -m cp pycistopic_output_wdl/DARs.pkl ~{output_dir}
        gsutil -m cp pycistopic_output_wdl/binarized_cell_topic.pkl ~{output_dir}
        gsutil -m cp pycistopic_output_wdl/binarized_topic_region.pkl ~{output_dir}
    >>>

    output {
        File DARs = 'pycistopic_output_wdl/DARs.pkl'
        File imputed_accessibility = 'pycistopic_output_wdl/Imputed_accessibility.pkl'
        File binarized_cell_topics = 'pycistopic_output_wdl/binarized_cell_topic.pkl'
        File binarized_topic_regions = 'pycistopic_output_wdl/binarized_topic_region.pkl'
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