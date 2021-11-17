#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "merges nearby DEL/DUP records within a certain window distance"

baseCommand: ["/bin/bash", "run_merge.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "apaul7/analysis:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "run_merge.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail
          INPUT="$1"
          OUTPUT="$2"
          DISTANCE="$3"
          /usr/local/bin/python3 /opt/git/merge-sv-records/merge.py -i "$INPUT" -o "$OUTPUT" -w "$DISTANCE"

          /usr/local/bin/bgzip "$OUTPUT"
          /usr/local/bin/tabix -p vcf "$OUTPUT".gz


inputs:
    input_vcf:
        type: File
        inputBinding:
            position: 1
    output_vcf_name:
        type: string?
        default: "record_merged.vcf"
        inputBinding:
            position: 2
    distance:
        type: int?
        default: 1000
        inputBinding:
            position: 3

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "$(inputs.output_vcf_name).gz"
        secondaryFiles: [.tbi]
