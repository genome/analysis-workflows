#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Concordance checking between Tumor and Normal BAM"
baseCommand: ["/usr/bin/somalier"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 8000
      tmpdirMin: 10000
    - class: DockerRequirement
      dockerPull: "brentp/somalier:v0.1.5"
arguments: ["-o", "concordance"]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-s"
            position: 1
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-f"
            position: 2
    cram_1:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [.crai]
    cram_2:
        type: File
        inputBinding:
            position: 4
        secondaryFiles: [.crai]
    cram_3:
        type: File?
        inputBinding:
            position: 5
        secondaryFiles: [.crai]
outputs:
    somalier_pairs:
        type: File
        outputBinding:
            glob: "concordance.somalier.pairs.tsv"
    somalier_samples:
        type: File
        outputBinding:
            glob: "concordance.somalier.samples.tsv"
