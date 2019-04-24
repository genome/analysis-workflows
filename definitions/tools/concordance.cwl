#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Concordance checking between Tumor and Normal BAM"
baseCommand: ["somalier"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 4000
      tmpdirMin: 10000
    - class: DockerRequirement
      dockerPull: "brentp/somalier"
arguments: ["-o", "concordance"]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-s"
            position: 1
    threads:
        type: string?
        inputBinding:
            prefix: "-t"
            position: 2
    min_depth:
        type: string?
        inputBinding:
            prefix: "-d"
            position: 3
    reference:
        type: string
        inputBinding:
            prefix: "-f"
            position: 4
    groups:
        type: File?
        inputBinding:
            prefix: "-g"
            position: 5
    ped:
        type: File?
        inputBinding:
            prefix: "-p"
            position: 6
    bam_1:
        type: File
        inputBinding:
            position: 7
        secondaryFiles: [.bai]
    bam_2:
        type: File
        inputBinding:
            position: 8
        secondaryFiles: [.bai]
outputs:
    somalier_pairs:
        type: File
        outputBinding:
            glob: "concordance.somalier.pairs.tsv"
    somalier_samples:
        type: File
        outputBinding:
            glob: "concordance.somalier.samples.tsv"
