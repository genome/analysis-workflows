#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "how_are_we_stranded_here"
baseCommand: ['check_strandedness']
arguments: [
    "--print_commands"
]

stdout: "$(inputs.reads1.nameroot).strandness_check.txt"

requirements:
    - class: ResourceRequirement
      ramMin: 16000
      tmpdirMin: 25000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "smk5g5/checkstranded:latest"
inputs:
    gtf_file:
        type: File
        inputBinding:
            prefix: --gtf
    kallisto_index:
        type: File
        inputBinding:
            prefix: --kallisto_index
    cdna_fasta:
        type: File
        inputBinding:
            prefix: --transcripts
    reads1:
        type: File
        inputBinding:
            prefix: "--reads_1"
    reads2:
        type: File
        inputBinding:
            prefix: "--reads_2"
outputs:
    check_strand:
        type: stdout
