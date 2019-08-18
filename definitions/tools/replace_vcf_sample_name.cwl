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

            if [[ "$1" =~ ".vcf.gz" ]];then
                zcat "$1" | sed -e "s/$2/$3/g" | /opt/htslib/bin/bgzip > sample_renamed.vcf.gz
            else
                cat "$1" | sed -e "s/$2/$3/g" | /opt/htslib/bin/bgzip > sample_renamed.vcf.gz
            fi

            /usr/bin/tabix -p vcf sample_renamed.vcf.gz

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
outputs:
    renamed_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputBinding:
            glob: sample_renamed.vcf.gz
