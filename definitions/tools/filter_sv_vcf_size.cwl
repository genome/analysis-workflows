#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "filter_sv_vcf_size.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.12"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "filter_sv_vcf_size.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail
          
          INPUT_VCF="$1"
          OUTPUT_VCF="$2"
          SV_SIZE="$3"
          FILTER_METHOD="$4"

          if [ "$FILTER_METHOD" == "max_len" ]; then
            echo "Running filter for max size svs"
            filter_expression="ABS(SVLEN) <= $SV_SIZE"
          elif [ "$FILTER_METHOD" ==  "min_len" ]; then
            echo "Running filter for min size svs"
            filter_expression="ABS(SVLEN) >= $SV_SIZE"
          else
            echo "Filter method: '$FILTER_METHOD' is not supported for size SV filtering"
            exit 1
          fi
          /opt/bcftools/bin/bcftools filter "$INPUT_VCF" -o "$OUTPUT_VCF" -i "$filter_expression"

inputs:
    input_vcf:
        type: File
        inputBinding:
            position: 1
        doc: "vcf file to filter"
    output_vcf_name:
        type: string?
        default: "filtered_sv_size.vcf"
        inputBinding:
            position: 2
        doc: "output vcf file name"
    sv_size:
        type: int?
        default: 50000
        inputBinding:
            position: 3
        doc: "size in bp used for filtering"
    size_method:
        type:
          - type: enum
            symbols: ["max_len", "min_len"]
        inputBinding:
          position: 4
        doc: "method for size filtering"

outputs:
    filtered_sv_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

