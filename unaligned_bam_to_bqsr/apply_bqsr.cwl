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
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/gatk-3.6:1"
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    reference:
        type: File
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
