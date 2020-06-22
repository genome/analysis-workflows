#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'apply BQSR'
baseCommand: ["/gatk/gatk", "--java-options", "-Xms14000m", "ApplyBQSR"]
arguments:
    ["-O", { valueFrom: $(runtime.outdir)/$(inputs.output_name).bam },
    "--static-quantized-quals", "10",
    "--static-quantized-quals", "20",
    "--static-quantized-quals", "30",
    "--create-output-bam-md5",
    "--use-original-qualities"]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.7.0"
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [.bai]
    bqsr_table:
        type: File
        inputBinding:
            prefix: "-bqsr"
            position: 3
    output_name:
        type: string?
        default: 'final'
outputs:
    bqsr_bam:
        type: File
        outputBinding:
            glob: $(inputs.output_name).bam
        secondaryFiles: [^.bai]