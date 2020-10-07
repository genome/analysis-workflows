#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "use the binomial/llr somatic filter to weed out low confidence variants"
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/somatic-llr-filter:v0.4.3
    - class: ResourceRequirement
      ramMin: 4000
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
    tumor_sample_name:
        type: string
        inputBinding:
            prefix: "--tumor-sample-name"
            position: -3
    normal_sample_name:
        type: string
        inputBinding:
            prefix: "--normal-sample-name"
    tumor_purity:
        type: float?
        inputBinding:
            prefix: "--tumor-purity"
            position: -4
        doc: "tumor cellularity fraction (range 0 to 1) - default 1"
    normal_contamination_rate:
        type: float?
        inputBinding:
            prefix: "--normal-contamination-rate"
            position: -5
        doc: "fraction of tumor present in the normal sample (range 0 to 1) - default 0"
outputs:
     somatic_llr_filtered_vcf:
         type: File
         outputBinding:
             glob: "somatic_llr_filtered.vcf"
