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
        
        ########### eRegulon correlation funcs ###########
        def p_adjust_bh(p: float):
            """
            Benjamini-Hochberg p-value correction for multiple hypothesis testing.
            from: pyCistopic: https://github.com/aertslab/pycisTopic/blob/d06246d9860157e028fcfee933fb3e784220b2c3/pycisTopic/diff_features.py#L747
            """
            p = np.asfarray(p)
            by_descend = p.argsort()[::-1]
            by_orig = by_descend.argsort()
            steps = float(len(p)) / np.arange(len(p), 0, -1)
            q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
            return q[by_orig]
        
        def eregulon_correlation(scplus_obj,
                         auc_key: str = 'eRegulon_AUC_filtered',
                         signature_key1: str = 'Gene_based',
                         signature_key2: str = 'Region_based',
                         nSignif: int = 3,
                         out_key: str = 'Unfiltered',
                         subset_cellids: List[str] = None
                        ):
            """
            Get correlation between gene-based and region-based eRegulon AUC
            Parameters
            ---------
            scplus_obj: :class:`SCENICPLUS`
                A :class:`SCENICPLUS` object
            auc_key: str, optional
                Name of the AUC matrix used to calculate the correlation, normally 'eRegulon_AUC'. 
                Must be a key in `scplus_obj.uns.keys()`.
            signature_key1: str, optional
                Variable used to calculate the correlation, normally 'Gene_based'. 
                Must be a key in `scplus_obj.uns[auc_key].keys()`.
            signature_key2: str, optional
                Variable used to calculate the correlation, normally 'Region_based'. 
                Must be a key in `scplus_obj.uns[auc_key].keys()`.
            nSignif: str, optional
                Number of digits to save.
            out_key : str, optional
                Ouput key. Correlations will be stored at `scplus_obj.uns['eRegulon_correlation'][out_key]`.
            subset_cellids: List, optional
                Subset of cells to be used to calculate correlations. Default: None (All)
            """

            gene_auc = scplus_obj.uns[auc_key][signature_key1].copy().T
            region_auc = scplus_obj.uns[auc_key][signature_key2].copy().T

            if subset_cellids is not None:
                cell_data = pd.DataFrame([x.rsplit('_', 1)[0] for x in gene_auc.columns],
                                        index=gene_auc.columns).iloc[:, 0]
                subset_cells = cell_data[cell_data.isin(subset_cellids)].index.tolist()
                gene_auc = gene_auc.loc[:, subset_cells]
                region_auc = region_auc.loc[:, subset_cells]

            # cistrome naming includes number of genes/regions, so need to create matching names
            # x.rsplit('_', 1) splits at first _ from the right
            gene_auc['id_short'] = gene_auc.index.map(lambda x: x.rsplit('_', 1)[0])
            gene_auc['id_full'] = gene_auc.index
            gene_auc = gene_auc.set_index('id_short')

            region_auc['id_short'] = region_auc.index.map(lambda x: x.rsplit('_', 1)[0])
            region_auc['id_full'] = region_auc.index
            region_auc = region_auc.set_index('id_short')

            corr_df = pd.DataFrame(columns=['TF', signature_key1, signature_key2, 'Rho', 'P-value'])

            for tf in gene_auc.index:
                # All TFs should match, but just in case
                if tf in region_auc.index:
                    # record orginal cistrome name for results
                    signature1_id = gene_auc.loc[tf, 'id_full']
                    signature2_id = region_auc.loc[tf, 'id_full']
                    # Exception in case TF is only expressed in 1 cell
                    # TFs expressed in few cells could be filtered too
                    try:
                        corr_1, _1 = st.pearsonr(gene_auc.loc[tf, gene_auc.columns != 'id_full'],
                                            region_auc.loc[tf, gene_auc.columns != 'id_full'])
                        x = {'TF': tf,
                            signature_key1: signature1_id,
                            signature_key2: signature2_id,
                            'Rho': round(corr_1,nSignif),
                            'P-value': _1}
                        corr_df = pd.concat([corr_df,
                                            pd.DataFrame(data=x, index=[0])],
                                            ignore_index=True)
                    except:
                        continue
            corr_df = corr_df.dropna()
            corr_df['Adjusted_p-value'] = p_adjust_bh(corr_df['P-value'])
            corr_df['Abs_rho'] = abs(corr_df['Rho'])
            corr_df.sort_values('Abs_rho', ascending=False, inplace=True)

            if not 'TF_cistrome_correlation' in scplus_obj.uns.keys():
                scplus_obj.uns['TF_cistrome_correlation'] = {}
            scplus_obj.uns['TF_cistrome_correlation'] = corr_df
        ############# end helper functions ###############

        # get eRegulon correlations
        eregulon_correlation(scplus_obj)

        thresholds = {
            'rho': [-0.75, 0.70],
            'n_targets': 0
        }

        # select for eRegulons 
        selected_cistromes = scplus_obj.uns['TF_cistrome_correlation'].loc[
            np.logical_or(
                scplus_obj.uns['TF_cistrome_correlation']['Rho'] > thresholds['rho'][1],
                scplus_obj.uns['TF_cistrome_correlation']['Rho'] < thresholds['rho'][0])]['Region_based'].to_list()

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