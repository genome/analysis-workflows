#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Replace the sample name in a gzipped vcf and return"
baseCommand: ['/bin/bash', 'replacer.sh']
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: ResourceRequirement
      ramMin: 8000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'replacer.sh'
        entry: |
            set -eou pipefail
            zcat $1 | sed -e "s/$2/$3/g" | /opt/htslib/bin/bgzip
arguments: [{valueFrom: $(inputs.input_vcf)}]
inputs:
    input_vcf:
        type: File
    sample_to_replace:
        type: string
        inputBinding:
            position: 2
    new_sample_name:
        type: string
        inputBinding:
            position: 3
stdout: "$(inputs.input_vcf.basename).renamed.gz"
outputs:
    renamed_vcf:
        type: stdout

