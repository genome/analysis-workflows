#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash","filter_sv_vcf_read_support.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "filter_sv_vcf_read_support.sh"
        entry: |
          set -eou pipefail
          
          input_vcf="$1"
          output_vcf="$2"
          paired_count="$3"
          paired_perc="$4"
          split_count="$5"
          split_perc="$6"
          vcf_source="$7"

          if [ "$vcf_source" == "smoove" ]; then
            echo "Running filter for smoove vcf"
            filter_expression="(AS >= $split_count && (AS / (RS+AS) >= $split_perc)) && (AP >= $paired_count && (AP / (AP+RP) >= $paired_perc))"
          elif [ "$vcf_source" ==  "manta" ]; then
            echo "Running filter for manta vcf"
            filter_expression="(SR[0:*]=\".\" || (SR[0:1] >= $split_count && (SR[0:1] / (SR[0:0]+SR[0:1]) >= $split_perc))) && (PR[0:1] >= $paired_count && (PR[0:1] / (PR[0:0]+PR[0:1]) >= $paired_perc))"
          else
            echo "vcf source: '$vcf_source' is not supported for SV filtering"
            exit 1
          fi

          /opt/bcftools/bin/bcftools filter "$input_vcf" -o "$output_vcf" -i "$filter_expression"

inputs:
    input_vcf:
        type: File
        inputBinding:
            position: 1
        doc: "vcf file to filter"
    output_vcf_name:
        type: string?
        default: "filtered_sv.vcf"
        inputBinding:
            position: 2
        doc: "output vcf file name"
    paired_count:
        type: int?
        default: 2
        inputBinding:
            position: 3
        doc: "number of alternate paired reads support needed to pass"
    paired_percentage:
        type: double?
        default: 0.2
        inputBinding:
            position: 4
        doc: "required alternate paired read abundance percentage to pass"
    split_count:
        type: int?
        default: 2
        inputBinding:
            position: 5
        doc: "if present in variant, number of alternate split reads support needed to pass"
    split_percentage:
        type: double?
        default: 0.2
        inputBinding:
            position: 6
        doc: "if present in variant, required alternate split read abundance percentage to pass"
    vcf_source:
        type:
          - type: enum
            symbols: ["manta", "smoove"]
        inputBinding:
          position: 7
        doc: "source caller of the vcf input file"

outputs:
    filtered_sv_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)
