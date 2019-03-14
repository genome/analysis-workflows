#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "filter vcf for variants with high percentage of mapq0 reads"
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/mapq0-filter:v0.1
    - class: ResourceRequirement
      ramMin: 8000
      tmpdirMin: 10000
arguments: 
    ["/bin/bash", "/usr/bin/mapq0_vcf_filter.sh",
    {valueFrom: "$(runtime.outdir)/mapq_filtered.vcf.gz"}]
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    tumor_bam:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [.bai]
    reference: 
         type: string
         inputBinding:
             position: 3
    threshold: 
         type: float
         inputBinding:
             position: 4
outputs:
     mapq0_filtered_vcf:
         type: File
         outputBinding:
             glob: "mapq_filtered.vcf.gz"
