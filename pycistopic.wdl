version 1.0

workflow pycisTopic {
    input {
    	String mallet_directory
        String cisTopic_obj_directory
        String tmp_directory
        String output_directory

        Array[Int] k_range = [3, 6, 9, 12]
        Int num_iter = 500
        Int seed = 555
        Int alpha = 50
        Boolean alpha_by_topic = true 
        Float eta = 0.1 
        Boolean eta_by_topic = false 
       
        #general parameters
        Int cpu = 8
        String memory = "64G"
        Int extra_disk_space = 0
        String docker = "mparikhbroad/cnmf:latest" # CHANGE
        Int preemptible = 0
    }

    # String output_directory_stripped = sub(output_directory, "/+$", "")

    call run_models {
            input:
                mallet_directory,
                cisTopic_obj_directory,
                tmp_directory,
                output_directory,
                
                k_range,
                num_iter,
                seed,
                alpha,
                alpha_by_topic,
                eta,
                eta_by_topic,

                cpu=cpu,
                memory=memory,
                extra_disk_space=extra_disk_space,
                docker=docker,
                preemptible=preemptible
        }
    
    output {
        File models = run_models.models
    }
}

task run_models {

    input {
        String mallet_directory
        String cisTopic_obj_directory
        String tmp_directory
        String output_directory

        Array[Int] k_range
        Int num_iter
        Int seed
        Int alpha
        Boolean alpha_by_topic
        Float eta 
        Boolean eta_by_topic
        
        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        source activate cnmf_env
        set -e

        mkdir -p outputs

        python <<CODE
        LDA_script = 'run_LDA.py'

        output_directory = 'outputs/'
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        run_name = 'LDA_run'

        prepare_cmd = ['python', cNMF_script, 'prepare', '--output-dir', output_directory, '--name', run_name, '-c', countfn, '-k', *k_range_array, '--n-iter', str(numiter), '--total-workers', str(numworkers), '--seed', str(seed), '--numgenes', str(numhvgenes), '--beta-loss', 'frobenius']
        print(' '.join(prepare_cmd), flush=True)
        check_call(prepare_cmd)
        CODE
        
        tar -zcf outputs.tar.gz outputs/cnmf_run
        
        tar -zcf cnmf.tar.gz outputs/cnmf_run
        rm -rf outputs/cnmf_run/cnmf_tmp
        gsutil -m rsync -r outputs/cnmf_run ~{output_dir}
    >>>

    output {
        File preparation_tar = 'outputs.tar.gz'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(adata_file, "GB")*4 + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}