#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash","filter_sv_vcf_read_support.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.12"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "filter_sv_vcf_read_support.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail

          abundance="$1"
          input_vcf="$2"
          output_vcf="$3"
          paired_count="$4"
          split_count="$5"
          vcf_source="$6"

          if [ "$vcf_source" == "smoove" ]; then
            echo "Running filter for smoove vcf"
            filter_expression="((AS >= $split_count) && (AP >= $paired_count) && ((AP+AS) / (AP+RP+AS+RS) >= $abundance))"
          elif [ "$vcf_source" ==  "manta" ]; then
            echo "Running filter for manta vcf"
            filter_expression="((SR[0:*]=\".\" || (SR[0:1] >= $split_count)) && (PR[0:1] >= $paired_count) && ((SR[0:*]=\".\" && (PR[0:1] / (PR[0:0]+PR[0:1]) >= $abundance))|| ((SR[0:1]+PR[0:1]) / (SR[0:0]+SR[0:1]+PR[0:1]+PR[0:0]) >= $abundance)))"
          else
            echo "vcf source: '$vcf_source' is not supported for SV filtering"
            exit 1
          fi

          /opt/bcftools/bin/bcftools filter "$input_vcf" -o "$output_vcf" -i "$filter_expression"

inputs:
    abundance_percentage:
        type: double?
        default: 0.1
        inputBinding:
            position: 1
        doc: "required alternate read abundance percentage to pass"
    input_vcf:
        type: File
        inputBinding:
            position: 2
        doc: "vcf file to filter"
    output_vcf_name:
        type: string?
        default: "filtered_sv.vcf"
        inputBinding:
            position: 3
        doc: "output vcf file name"
    paired_count:
        type: int?
        default: 2
        inputBinding:
            position: 4
        doc: "number of alternate paired reads support needed to pass"
    split_count:
        type: int?
        default: 2
        inputBinding:
            position: 5
        doc: "if present in variant, number of alternate split reads support needed to pass"
    vcf_source:
        type:
          - type: enum
            symbols: ["manta", "smoove"]
        inputBinding:
          position: 6
        doc: "source caller of the vcf input file"

outputs:
    filtered_sv_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)
