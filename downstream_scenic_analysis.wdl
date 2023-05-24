version 1.0

workflow filter_scenicplus {
    call run_filter_scplus

    output {
        File downstream_analysis_output = run_filter_scplus.filtered_scplus_object
    }
}

task run_filter_scplus {
    input {
        String output_dir # gbucket (add / at end)
        File scplus_obj_file # scplus_obj2.pkl
        File region_ranking_file # region_ranking.pkl
        File gene_ranking_file # gene_ranking.pkl

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:1.0.0"
        Int preemptible = 0
        Int disk_space = 128
    }

    command <<<
        set -e 
        
        mkdir tmpdir
        mkdir downstream_scenicplus_output_wdl

        python << CODE
        # imports
        import os
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import pickle
        import dill
        from pycisTopic.cistopic_class import *
        import scipy.stats as st

        # load scenic plus object
        scplus_obj = dill.load(open('~{scplus_obj_file}', 'rb'))

        # filter scenicplus output
        from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
        apply_std_filtering_to_eRegulons(scplus_obj)    

        # score the enrichment of eRegulons using the AUCell function
        from scenicplus.eregulon_enrichment import score_eRegulons
        region_ranking = dill.load(open('~{region_ranking_file}', 'rb')) 
        gene_ranking = dill.load(open('~{gene_ranking_file}', 'rb')) 
        score_eRegulons(scplus_obj,
                        ranking = region_ranking,
                        eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                        key_added = 'eRegulon_AUC_filtered',
                        enrichment_type= 'region',
                        auc_threshold = 0.05,
                        normalize = False,
                        n_cpu = 5)
        score_eRegulons(scplus_obj,
                        gene_ranking,
                        eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                        key_added = 'eRegulon_AUC_filtered',
                        enrichment_type = 'gene',
                        auc_threshold = 0.05,
                        normalize= False,
                        n_cpu = 5)
        
        # eRegulon dimensionality reduction
        from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
        run_eRegulons_umap(
            scplus_obj = scplus_obj,
            auc_key = 'eRegulon_AUC',
            reduction_name = 'eRegulons_UMAP', #overwrite previously calculated UMAP
        )
        run_eRegulons_tsne(
            scplus_obj = scplus_obj,
            auc_key = 'eRegulon_AUC',
            reduction_name = 'eRegulons_tSNE', #overwrite previously calculated tSNE
        )
        
        # compute correlation between TF expression and target region enrichment scores (AUC values)
        from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

        generate_pseudobulks(
                scplus_obj = scplus_obj,
                variable = 'GEX_states',
                auc_key = 'eRegulon_AUC_filtered',
                signature_key = 'Gene_based')
        generate_pseudobulks(
                scplus_obj = scplus_obj,
                variable = 'GEX_states',
                auc_key = 'eRegulon_AUC_filtered',
                signature_key = 'Region_based')

        TF_cistrome_correlation(
                    scplus_obj,
                    use_pseudobulk = True,
                    variable = 'GEX_states',
                    auc_key = 'eRegulon_AUC_filtered',
                    signature_key = 'Gene_based',
                    out_key = 'filtered_gene_based')
        TF_cistrome_correlation(
                    scplus_obj,
                    use_pseudobulk = True,
                    variable = 'GEX_states',
                    auc_key = 'eRegulon_AUC_filtered',
                    signature_key = 'Region_based',
                    out_key = 'filtered_region_based')

        # select for eRegulons 
        thresholds = {
                'rho': [-0.65, 0.60],
                'n_targets': 0
        }
        
        selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
                np.logical_or(
                        scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
                        scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
                )]['Cistrome'].to_list()
        selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
        selected_eRegulons_gene_sig = [
                x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
                if x.split('_(')[0] in selected_eRegulons]
        selected_eRegulons_region_sig = [
                x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
                if x.split('_(')[0] in selected_eRegulons]


        # save the results in the scenicplus object
        scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig, 'Region_based': selected_eRegulons_region_sig}
        
        # save new scplus object
        dill.dump(scplus_obj, open(os.path.join('downstream_scenicplus_output_wdl', 'filtered_scplus_obj.pkl'), 'wb'), protocol=-1)
                
        CODE

        gsutil -m cp downstream_scenicplus_output_wdl/filtered_scplus_obj.pkl ~{output_dir}
    >>>

    output {
        File filtered_scplus_object = 'downstream_scenicplus_output_wdl/filtered_scplus_obj.pkl'
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