#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Sanitize mutect vcfs, which have format tags containing vcf-standard illegal underscores"
baseCommand: ['/bin/bash', 'sanitizer.sh']
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: ResourceRequirement
      ramMin: 8000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'sanitizer.sh'
        entry: |
            set -eou pipefail
            zcat $1 | sed -e "s/ALT_F1R2/ALTF1R2/g" | sed -e "s/ALT_F2R1/ALTF2R1/g" | sed -e "s/REF_F1R2/REFF1R2/g" | sed -e "s/REF_F2R1/REFF2R1/g" | /opt/htslib/bin/bgzip
arguments: [{valueFrom: $(inputs.input_vcf)}]
inputs:
    input_vcf:
        type: File
stdout: "$(inputs.input_vcf.basename).sanitized.gz"
outputs:
    sanitized_vcf:
        type: stdout

