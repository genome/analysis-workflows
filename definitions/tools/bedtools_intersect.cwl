#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/usr/local/bin/bedtools", "intersect"]

requirements:
    - class: DockerRequirement
      dockerPull: "ernfrid/bedtools:v2.27.1"
    - class: ResourceRequirement
      ramMin: 4000

stdout: "$(inputs.output_name)"

inputs:
    file_a:
        type: File
        inputBinding:
            position: 1
            prefix: "-a"
        doc: "BAM/BED/GFF/VCF file to compare to file b"
    file_b:
        type: File
        inputBinding:
            position: 2
            prefix: "-b"
        doc: "BAM/BED/GFF/VCF file to compare to file a"
    output_file_a:
        type: boolean?
        default: true
        inputBinding:
            position: 3
            prefix: "-wa"
        doc: "Write the original entry in A for each overlap, default to true"
    output_header:
        type: boolean?
        default: true
        inputBinding:
            position: 4
            prefix: "-header"
        doc: "Print the header from the A file prior to results, default to true"
    output_name:
        type: string
    unique_result:
        type: boolean?
        default: true
        inputBinding:
            position: 5
            prefix: "-u"
        doc: "Write original A entry once if any overlaps found in B, default to true"

outputs:
    intersect_result:
        type: stdout
