#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard MergeVcfs"
baseCommand: ["/usr/bin/java", "-jar", "/opt/picard/picard.jar", "MergeVcfs"]
requirements:
    - class: ResourceRequirement
      ramMin: 40000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
arguments:
    ["O=merged.vcf.gz"]
inputs:
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
            glob: "merged.vcf.gz"
        secondaryFiles: [.tbi]
