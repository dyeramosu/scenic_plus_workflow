version 1.0

workflow cisTarget {
    call run_pycistarget

    output {
        File pycistarget_output = run_pycistarget.pycistarget_object
    }
}

task run_pycistarget {
    input {
        String output_dir # gbucket
        File binarized_topic_region_file
        File DARs_file
        File rankings_db_file # hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
        File scores_db_file # hg38_screen_v10_clust.regions_vs_motifs.scores.feather
        File motif_annotation_file # motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

        Int cpu = 24
        Int memory = 256
        String docker = "dyeramosu/scenic_plus_terra:latest"
        Int preemptible = 0
        Int disk_space = 128
    }

    command <<<
        set -e 
        
        mkdir tmpdir
        mkdir pycistarget_output_wdl

        python << CODE
        # imports
        import os
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import scanpy as sc
        import pickle
        from pycisTopic.cistopic_class import *

        # load candidate enhancer regions
        region_bin_topics = pickle.load('~{binarized_topic_region_file}', 'rb')
        markers_dict = pickle.load('~{DARs_file}', 'rb')
        
        markers_dict.pop('Unknown')

        # convert to dictionary of pyranges objects
        import pyranges as pr
        from pycistarget.utils import region_names_to_coordinates

        region_sets = {}
        region_sets['topics'] = {}
        region_sets['DARs'] = {}

        for topic in region_bin_topics.keys():
            # only keep regions on known chromosomes
            regions = region_bin_topics[topic].index[region_bin_topics[topic].index.str.startswith('chr')] 
            region_sets['topics'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

        for DAR in markers_dict.keys():
            # only keep regions on known chromosomes
            regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] 
            region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

        # define rankings, score and motif annotation database
        rankings_db = '~{rankings_db_file}'
        scores_db = '~{scores_db_file}'
        motif_annotation = '~{motif_annotation_file}'

        # run pycistarget using the run_pycistarget wrapper function
        from scenicplus.wrappers.run_pycistarget import run_pycistarget

        run_pycistarget(
            region_sets = region_sets,
            species = 'homo_sapiens',
            save_path = 'pycistarget_output_wdl/',
            ctx_db_path = rankings_db,
            dem_db_path = scores_db,
            path_to_motif_annotations = motif_annotation,
            run_without_promoters = True,
            n_cpu = 1,
            _temp_dir = 'tmpdir/',
            annotation_version = 'v10nr_clust',
            )
        
        # try to open menr file
        import dill
        menr = dill.load(open(os.path.join('pycistarget_output_wdl', 'menr.pkl'), 'rb'))
       
        CODE

        tar -czvf pycistarget_output.tar.gz pycistarget_output_wdl
        gsutil rsync -r pycistarget_output_wdl ~{output_dir}
    >>>

    output {
        File pycistarget_object = 'pycistarget_output.tar.gz'
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