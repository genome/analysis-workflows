#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit vcf2bed wrapper"
baseCommand: ["/usr/bin/bsvcf2bed"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 2
#Creates a gzipped bed and a bedgraph that leaves out MT, random, GL contigs, etc
arguments: [
    "$(runtime.outdir)/cpgs.bed.gz",
    "$(runtime.outdir)/cpgs.bedgraph"
]
inputs:
    vcf:
        type: File
        inputBinding:
            position: -1
outputs:
    cpgs:
        type: File
        outputBinding:
            glob: "cpgs.bed.gz"
    cpg_bedgraph:
        type: File
        outputBinding:
            glob: "cpgs.bedgraph"
