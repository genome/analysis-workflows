#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Subset gvcf on Intervals File"
baseCommand: ["/gatk/gatk", "-jar", "SelectVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 20000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.7.0"
arguments:
    ["--remove-unused-alternates", "--preserve-alleles", "-O interval_subset.vcf.gz"]
inputs:
    merged_gvcf:
        type: File
        inputBinding:
            prefix: "-V"
            position: 0
    interval_list:
        type: File
        inputBinding: 
            prefix: "-L"
            position: 1
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 2
    
outputs:
    interval_gvcf:
        type: File
        outputBinding:
            glob: "interval_subset.vcf.gz"
        secondaryFiles: [.tbi]