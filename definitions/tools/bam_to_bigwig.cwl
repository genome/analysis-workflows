#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'cgpBigWig Converting BAM to BigWig'
baseCommand: ["bam2bw", "-F", "1024"]
arguments: ["-o", { valueFrom: $(inputs.bam.nameroot).bw}]
requirements:
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/cgpbigwig:1.4.0--h93d22ca_0"
    - class: ResourceRequirement
      ramMin: 32000
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
            prefix: '-i'
        secondaryFiles: [^.bai, .bai]
    reference:
        type:
            - string
            - File
        inputBinding:
            position: 2
            prefix: '-r'
        secondaryFiles: [.fai, ^.dict]

outputs:
    outfile:
        type: File
        outputBinding:
            glob: $(inputs.bam.nameroot).bw
