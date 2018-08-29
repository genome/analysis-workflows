#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "filter_vcf for variants with high percentage of mapq0 reads"
arguments: [
    "/bin/bash /usr/bin/mapq0_vcf_filter.sh",
    {valueFrom: "$(runtime.outdir)/mapq_filtered.vcf"}]
]
inputs:
    vcf:
        type: File
        inputBinding:
            position: -4
     tumor_bam: 
        type: File
        inputBinding:
            position: -3
    reference: 
        type: string
        inputBinding:
            position: -2
    threshold: 
        type: float
        inputBinding:
            position: -1

outputs:
    mapq0_filtered_vcf:
        type: File
        outputBinding:
            glob: "mapq_filtered.vcf"
