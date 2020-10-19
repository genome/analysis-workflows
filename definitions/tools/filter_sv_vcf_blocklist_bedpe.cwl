#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "filter_sv_vcf_blocklist_bedpe.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "mgibio/basespace_chromoseq:v12"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "filter_sv_vcf_blocklist_bedpe.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail
          set -o errexit
          
          INPUT_VCF="$1"
          BL_BEDPE="$2"
          SLOPE="$3"
          OUT_BASE="$4"
          
          #BASE=`basename $1`
          #NAMEROOT=`echo $BASE | perl -pe 's/.vcf(.gz)?$//g'`
          
          if [[ "$BL_BEDPE" == 'NONE' ]]; then
              /usr/local/bin/bedtools sort -header -i "$INPUT_VCF" > $OUT_BASE.vcf 
              /opt/htslib/bin/bgzip $OUT_BASE.vcf
              /usr/bin/tabix -p vcf $OUT_BASE.vcf.gz
          else
              #CNVkit outputs invalid format like CIPOS=.,894;CIEND=.,894, which can cause svtools vcftobedpe fail
              if [[ "$INPUT_VCF" =~ \.vcf\.gz$ ]]; then
                  /bin/zcat "$INPUT_VCF" | /bin/sed -E 's/CIPOS=\.,[0-9]+;CIEND=\.,[0-9]+/CIPOS=0,0;CIEND=0,0/g' > fixed_input.vcf
              else
                  /bin/sed -E 's/CIPOS=\.,[0-9]+;CIEND=\.,[0-9]+/CIPOS=0,0;CIEND=0,0/g' "$INPUT_VCF" > fixed_input.vcf
              fi
              #svtools vcftobedpe can take either .vcf or .vcf.gz
              /opt/conda/envs/python2/bin/svtools vcftobedpe -i fixed_input.vcf -o tmp.bedpe
              /bin/grep '^#' tmp.bedpe > tmp.header
              /usr/local/bin/bedtools pairtopair -is -slop "$SLOPE" -type notboth -a tmp.bedpe -b "$BL_BEDPE" | /bin/cat tmp.header /dev/stdin | /opt/conda/envs/python2/bin/svtools bedpetovcf -i /dev/stdin | /opt/conda/envs/python2/bin/svtools vcfsort /dev/stdin > "$OUT_BASE.vcf"
          
              /opt/htslib/bin/bgzip $OUT_BASE.vcf
              /usr/bin/tabix -p vcf $OUT_BASE.vcf.gz
          fi

inputs:
    input_vcf:
        type: File
        inputBinding:
            position: 1
        doc: "vcf file to filter"
    blocklist_bedpe:
        type: string
        inputBinding:
            position: 2
        doc: "blocklist bedpe file"
    slope:
        type: int?
        default: 100
        inputBinding:
            position: 3
        doc: "slope(bp) used in bedtools pairtopair"
    output_vcf_basename:
        type: string?
        default: "blocklist_filtered"
        inputBinding:
            position: 4
        doc: "output vcf basename"
outputs:
    filtered_sv_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputBinding:
            glob: $(inputs.output_vcf_basename).vcf.gz
