#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "filter_sv_vcf_depth.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "filter_sv_vcf_depth.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail
          
          INPUT_VCF="$1"
          OUTPUT_VCF="$2"
          DEL_DEPTH="$3"
          DUP_DEPTH="$4"
          VCF_SOURCE="$5"

          if [ "$VCF_SOURCE" == "cnvnator" ]; then
            echo "Running filter for cnvnator vcf"
            filter_expression="(natorRD>$DUP_DEPTH | natorRD<$DEL_DEPTH)"
          elif [ "$VCF_SOURCE" ==  "cnvkit" ]; then
            echo "Running filter for cnvkit vcf"
            filter_expression="(FOLD_CHANGE>$DUP_DEPTH | FOLD_CHANGE<$DEL_DEPTH)"
          else
            echo "vcf source: '$VCF_SOURCE' is not supported for SV filtering"
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
        default: "filtered_cnv.vcf"
        inputBinding:
            position: 2
        doc: "output vcf file name"
    deletion_depth:
        type: double?
        default: 0.75
        inputBinding:
            position: 3
        doc: "Depth for deletions must be less than this cutoff to pass the filter."
    duplication_depth:
        type: double?
        default: 1.25
        inputBinding:
            position: 4
        doc: "Cutoff depth for duplications, greater than passes filter"
    vcf_source:
        type:
          - type: enum
            symbols: ["cnvnator", "cnvkit"]
        inputBinding:
          position: 5
        doc: "source caller of the vcf input file"

outputs:
    filtered_sv_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)
