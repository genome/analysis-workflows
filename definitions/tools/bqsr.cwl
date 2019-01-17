#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'create BQSR table'
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator"]
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/bqsr.table },
    "--preserve_qscores_less_than", "6",
    "--disable_auto_index_creation_and_locking_when_reading_rods",
    "--disable_bam_indexing",
    "-dfrac", ".1",
    "-nct", "4"]
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
            position: 2
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 3
        secondaryFiles: [.bai]
    known_sites:
        type:
            type: array
            items: File
            inputBinding:
                prefix: "-knownSites"
                position: 4
        secondaryFiles: [.tbi]
    intervals:
        type:
            type: array
            items: string
            inputBinding:
                prefix: "-L"
        inputBinding:
            position: 1
        default: [chr1, chr2, chr3, chr4, chr5,chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22]
outputs:
    bqsr_table:
        type: File
        outputBinding:
            glob: "bqsr.table"
