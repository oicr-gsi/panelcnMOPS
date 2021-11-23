version 1.0

workflow PanelCN {
    meta {
        author: "Vivek Alamuri"
        email: "valamuri@oicr.on.ca"
        description: "Workflow that runs Cn mops R package"
        dependencies: [
            {
                name: "rstats/4.0",
                url: "https://www.r-project.org/"
            }
        ]
    }

    input {
        Array[File] ref_bams_files
        Array[File] test_bams_files
        File bed_file
        String output_directory
        File script
    }

    call runPanelCN {
        input:
            ref_bams_files = ref_bams_files,
            test_bams_files = test_bams_files,
            bed_file = bed_file,
            output_directory = output_directory,
            script = script
    }
}

task runPanelCN {
    input {
        Array[File] ref_bams_files
        Array[File] test_bams_files
        File bed_file
        String output_directory
        File script

        String modules = "rstats/4.0"
        Int threads = 8
        Int jobMemory = 64
        Int timeout = 72
    }

    parameter_meta {
        ref_bams_files: ""
        test_bams_files: ""
        bed_file: ""
        output_directory: "output directory"
        script: "path to R script"
        modules: "Names and versions of modules to load"
        threads: "Requested CPU threads"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }

    command <<< 
        Rscript ~{script} -t ~{test_bams_files} -r {ref_bams_files} -b ~{bed_file} \
        -o ~{output_directory}
    >>>
    
    runtime {
        memory:  "~{jobMemory} GB"
        modules: "~{modules}"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }

}