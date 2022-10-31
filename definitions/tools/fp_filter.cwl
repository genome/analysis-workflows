#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "False Positive filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/fpfilter.pl"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/fp_filter-cwl:1.0.1"
    - class: InitialWorkDirRequirement
      listing:
      - $(inputs.ref_fai)
      - $(inputs.ref_dict)
      - $(inputs.reference)
arguments:
    ["--bam-readcount", "/usr/bin/bam-readcount",
    "--samtools", "/opt/samtools/bin/samtools",
    "--output", { valueFrom: $(runtime.outdir)/$(inputs.output_vcf_basename).vcf }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "--reference"
            position: 1
    ref_fai:
        type: File
    ref_dict:
        type: File
    bam:
        type: File
        inputBinding:
            prefix: "--bam-file"
            position: 2
    vcf:
        type: File
        inputBinding:
            prefix: "--vcf-file"
            position: 3
    output_vcf_basename:
        type: string?
        default: fpfilter
    sample_name:
        type: string?
        default: 'TUMOR'
        inputBinding:
            prefix: "--sample"
            position: 4
    min_var_freq:
        type: float?
        default: 0.05
        inputBinding:
            prefix: "--min-var-freq"
            position: 5
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_basename).vcf

