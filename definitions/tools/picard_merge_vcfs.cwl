#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard MergeVcfs"
baseCommand: ["/usr/bin/java", "-jar", "/opt/picard/picard.jar", "MergeVcfs"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
arguments:
    ["O=$(inputs.merged_vcf_basename).vcf.gz"]
inputs:
    merged_vcf_basename:
        type: string?
        default: "merged"
    gvcfs:
        type:
            type: array
            items: File
            inputBinding:
                prefix: "I="
                separate: false
        inputBinding:
            position: 0
outputs:
    merged_vcf:
        type: File
        outputBinding:
            glob: $(inputs.merged_vcf_basename).vcf.gz
        secondaryFiles: [.tbi]
