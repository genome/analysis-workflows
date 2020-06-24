#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Knight lab bamMetrics"
baseCommand: ["/gscmnt/gc2698/jin810/bamMetrics_docker"]
requirements:
    - class: ResourceRequirement
      ramMin: 10000
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
      
inputs:
    bed:
        type: File
        inputBinding:
            prefix: "-b"
            position: 1
    reference:
        type:
            - string 
            - File
        inputBinding:
            prefix: "-r"
            position: 2
    output_filename:
       type: string
       default: "bam_metrics.txt"
       inputBinding:
           prefix: "-o"
           position: 3
    cram:
        type: File
        inputBinding:
            position: 4
outputs:
    bam_metrics:
        type: File
        outputBinding:
            glob: $(inputs.output_filename)