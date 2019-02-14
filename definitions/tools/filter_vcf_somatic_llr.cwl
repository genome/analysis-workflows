#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "use the binomial/llr somatic filter to weed out low confidence variants"
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/somatic-llr-filter:v0.1
baseCommand: ["/opt/conda/bin/python3","/usr/bin/somatic_llr_filter.py"]
arguments: 
    ["--overwrite", #we expect to have to overwrite the SOMATIC field
    {valueFrom: "$(runtime.outdir)/somatic_llr_filtered.vcf"}]
inputs:
    vcf:
        type: File
        inputBinding:
            position: -1
    threshold: 
         type: float
         inputBinding:
             prefix: "--llr-threshold"
             position: -2
outputs:
     somatic_llr_filtered_vcf:
         type: File
         outputBinding:
             glob: "somatic_llr_filtered.vcf"
