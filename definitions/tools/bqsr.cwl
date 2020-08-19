#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'create BQSR table'
baseCommand: ["/gatk/gatk","--java-options", "-Xmx16g", "BaseRecalibrator"]
arguments:
    ["-O", { valueFrom: $(runtime.outdir)/bqsr.table }
    ]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
                prefix: "--known-sites"
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
