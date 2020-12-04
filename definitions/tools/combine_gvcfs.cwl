#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK CombineGVCFs"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx8g", "CombineGVCFs"]
requirements:
    - class: ResourceRequirement
      ramMin: 10000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "--reference"
            position: 1
    gvcfs:
        type:
          type: array
          items: File
          inputBinding:
            prefix: "--variant"
        inputBinding:
            position: 2
    output_file_name:
        type: string
        default: "merged.g.vcf.gz"
        inputBinding:
            prefix: "--output"
            position: 3
outputs:
    gvcf:
        type: File
        outputBinding:
            glob: $(inputs.output_file_name)
        secondaryFiles: [.tbi]
