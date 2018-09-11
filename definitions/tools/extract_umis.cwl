#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'extract umis from bam'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/fgbio-0.5.0.jar", "ExtractUmisFromBam"]
arguments:
    ["--molecular-index-tags", "ZA", "ZB", "--single-tag", "RX",
    "--output", { valueFrom: "$(runtime.outdir)/umi_extracted.bam"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
    read_structure:
        type: string[]
        inputBinding:
            prefix: "--read-structure"
outputs:
    umi_extracted_bam:
        type: File
        outputBinding:
            glob: "umi_extracted.bam"
