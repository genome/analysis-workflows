#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "convert a biscuit pileup vcf to a bed and bedgraph file, via a wrapper script"
baseCommand: ["/usr/bin/bsvcf2bed"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 2
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite:v1.3"
#Creates a gzipped bed and a bedgraph that leaves out MT, random, GL contigs, etc
arguments: [
    "$(runtime.outdir)/cpgs.bed.gz",
    "$(runtime.outdir)/cpgs.bedgraph"
]
inputs:
    vcf:
        type: File
        inputBinding:
            position: -2
    reference:
        type: string
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
