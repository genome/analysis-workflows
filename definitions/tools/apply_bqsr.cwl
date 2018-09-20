#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'apply BQSR'
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "PrintReads"]
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/Final.bam },
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
        type: string
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
            prefix: "-BQSR"
            position: 3
outputs:
    bqsr_bam:
        type: File
        outputBinding:
            glob: "Final.bam"
        secondaryFiles: [^.bai]
