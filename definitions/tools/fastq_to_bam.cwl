#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'fastq to bam'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/picard/picard.jar", "FastqToSam"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/unaligned.bam }]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment
inputs:
    read1_fastq:
        type: File
        inputBinding:
            prefix: "FASTQ="
            position: 1
    read2_fastq:
        type: File
        inputBinding:
            prefix: "FASTQ2="
            position: 2
    sample_name:
        type: string
        inputBinding:
            prefix: "SAMPLE_NAME="
            position: 3
    library_name:
        type: string
        inputBinding:
            prefix: "LIBRARY_NAME="
            position: 4
    platform_unit:
        type: string
        inputBinding:
            prefix: "PLATFORM_UNIT="
            position: 5
    platform:
        type: string
        inputBinding:
            prefix: "PLATFORM="
            position: 6
outputs:
    bam:
        type: File
        outputBinding:
            glob: "unaligned.bam"
