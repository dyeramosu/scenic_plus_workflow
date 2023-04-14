version 1.0

workflow cisTarget {
    call run_pycistarget

    output {

    }
}

task run_pycistarget {
    input {

    }

    command {
        set -e 
        mkdir tmp
        mkdir pycistarget_output_wdl

        

    }

    output {

    }

    runtime {

    }
}