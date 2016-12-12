#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bwa_mem"
baseCommand: ["/usr/local/bin/bwa", "mem"]
arguments: 
    ["-K", "10000000",
    "-t", "8",
    "-Y"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/bwa-0.7.15:1"
stdout: 'refAlign.bam'
inputs:
    readgroup:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    reference:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac]
    fastq:
        type: File
        inputBinding:
            position: 3
    fastq2:
        type: File
        inputBinding:
            position: 4
outputs:
    aligned_bam:
        type: stdout
