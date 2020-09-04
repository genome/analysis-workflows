#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'apply BQSR'
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "PrintReads"]
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/$(inputs.output_name).cram },
    "-preserveQ", "6",
    "-SQQ", "10",
    "-SQQ", "20",
    "-SQQ", "30",
    "-nct", "8",
    "--disable_indel_quals"]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.6.0"
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    cram:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [.crai]
    bqsr_table:
        type: File
        inputBinding:
            prefix: "-BQSR"
            position: 3
    output_name:
        type: string?
        default: 'final'
outputs:
    bqsr_cram:
        type: File
        outputBinding:
            glob: $(inputs.output_name).cram
        secondaryFiles: [^.crai]
