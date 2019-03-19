#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard: Revert Input"
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/picard/picard.jar", "RevertSam", "VALIDATION_STRINGENCY=SILENT"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq"
arguments:
    - valueFrom: "reverted.bam"
      position: 2
      prefix: "O="

inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
            separate: false
            position: 1
outputs:
    reverted_bam:
        type: File
        outputBinding:
            glob: "reverted.bam"
