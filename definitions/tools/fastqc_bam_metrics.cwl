#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: "v1.0"
label: "Report metrics on bam files using fastqc"
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 2
    - class: DockerRequirement
      dockerPull: mgibio/fastqc:0.11.9

baseCommand: ["/usr/bin/perl", "/usr/local/bin/FastQC/fastqc"]
arguments: ["-q", "--extract", "-o", "$(runtime.outdir)"]

inputs:
  bam:
    type: File
    inputBinding:
        position: 1

outputs:
    fastqc_all_data:
      type: Directory
      outputBinding:
        glob: "$(inputs.bam.nameroot)_fastqc"

    fastqc_txt_report:
      type: File
      outputBinding:
        glob: "$(inputs.bam.nameroot)_fastqc/fastqc_data.txt"
