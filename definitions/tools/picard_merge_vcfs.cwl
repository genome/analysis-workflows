#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard MergeVcfs"
baseCommand: ["/usr/bin/java", "-jar", "-Xmx38g", "/opt/picard/picard.jar", "MergeVcfs"]
requirements:
    - class: ResourceRequirement
      ramMin: 40000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "broadinstitute/picard:2.23.6"
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
