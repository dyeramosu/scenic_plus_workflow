version 1.0

workflow get_scores {
    call run_get_gene_activity

    output {
        File get_scores_output = run_get_gene_activity.scores_object
    }
}

task run_get_gene_activity {
    input {
        String output_dir # gbucket (no / at end)
        File DARs_file # DARs.pkl
        File cistopic_file # cistopic_obj_plus_model.pkl
        File imputed_acc_file # Imputed_accessibility.pkl 

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:1.0.0"
        Int preemptible = 0
        Int disk_space = 128
    }

    command <<<
        set -e 
        
        mkdir tmpdir
        mkdir get_scores_output_wdl

        python << CODE
        # imports
        import os
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import pickle
        from pycisTopic.cistopic_class import *

        # LOAD OBJECTS
        cistopic_obj = pickle.load(open('~{cistopic_file}', 'rb'))
        imputed_acc_obj = pickle.load(open('~{imputed_acc_file}', 'rb'))
        DARs_dict = pickle.load(open('~{DARs_file}', 'rb'))

        # RETRIEVE GENE ANNOTATION AND CHROMOSOME SIZES FOR OUR GENOME
        ## get TSS annotations
        import pybiomart as pbm
        import pyranges as pr

        dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://sep2019.archive.ensembl.org/')
        annot = dataset.query(attributes=['chromosome_name', 'start_position', 'end_position', 'strand', 'external_gene_name', 
                                        'transcription_start_site', 'transcript_biotype'])
        annot['Chromosome/scaffold name'] = 'chr' + annot['Chromosome/scaffold name'].astype(str)
        annot.columns=['Chromosome', 'Start', 'End', 'Strand', 'Gene','Transcription_Start_Site', 'Transcript_type']
        annot = annot[annot.Transcript_type == 'protein_coding']
        annot.Strand[annot.Strand == 1] = '+'
        annot.Strand[annot.Strand == -1] = '-'
        pr_annotation = pr.PyRanges(annot.dropna(axis = 0))

        ## get chromosome sizes
        import pandas as pd
        import requests
        target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
        chromsizes=pd.read_csv(target_url, sep='\t', header=None)
        chromsizes.columns=['Chromosome', 'End']
        chromsizes['Start']=[0]*chromsizes.shape[0]
        chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
        chromsizes=pr.PyRanges(chromsizes)

        # INFER GENE ACTIVITY
        from pycisTopic.gene_activity import *
        gene_act, weights = get_gene_activity(imputed_acc_obj, # Region-cell probabilities
                        pr_annotation, # Gene annotation
                        chromsizes, # Chromosome size
                        use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
                        upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                                            #these bp will be taken (1kbp here)
                        downstream=[1000,100000], # Search space downstream
                        distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
                        decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
                        extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                                            #this weight)
                        extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
                        gene_size_weight=False, # Whether to add a weights based on the length of the gene
                        gene_size_scale_factor='median', # Dividend to calculate the gene size weight. Default is the median value of all genes
                                            #in the genome
                        remove_promoters=False, # Whether to remove promoters when computing gene activity scores
                        average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                                            #activity score
                        scale_factor=1, # Value to multiply for the final gene activity matrix
                        extend_tss=[10,10], # Space to consider a promoter
                        gini_weight = True, # Whether to add a gini index weight. The more unique the region is, the higher this weight will be
                        return_weights= True, # Whether to return the final weights
                        project='Gene_activity') # Project name for the gene activity object 
        
        # INFER DAGs 
        markers_dict = find_diff_features(cistopic_obj,
                      gene_act,
                      variable='states',
                      var_features=None,
                      contrasts=None,
                      adjpval_thr=0.05,
                      log2fc_thr=np.log2(1.5),
                      n_cpu=5,
                      _temp_dir=os.path.abspath("tmpdir/"),
                      split_pattern = '-')
        
        # SAVE
        with open(os.path.join('get_scores_output_wdl', 'Gene_activity.pkl'), 'wb') as f:
            pickle.dump(gene_act, f)
        
        with open(os.path.join('get_scores_output_wdl', 'weights.pkl'), 'wb') as f:
            pickle.dump(weights, f)

        with open(os.path.join('get_scores_output_wdl', 'DAGs.pkl'), 'wb') as f:
            pickle.dump(markers_dict, f)
       
        CODE

        tar -czvf get_scores_output.tar.gz get_scores_output_wdl
        gsutil rsync -r get_scores_output_wdl ~{output_dir}
    >>>

    output {
        File scores_object = 'get_scores_output.tar.gz'
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