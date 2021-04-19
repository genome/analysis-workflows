#!/usr/bin/env cwl-runner
 
cwlVersion: v1.0
class: CommandLineTool
label: "Adds an INFO tag (VALIDATED) flagging variants in the pipeline vcf present in a previously validated vcf file"

requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: ResourceRequirement
      ramMin: 8000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'annotate.sh'
        entry: |
            set -eou pipefail

            PIPELINE_VCF="$1"

            if [ "$#" -eq 2 ]; then
                VALIDATED_VCF="$2"
                /opt/bcftools/bin/bcftools view -f PASS -Oz -o pass_filtered_validated_variants.vcf.gz $VALIDATED_VCF
                /opt/bcftools/bin/bcftools index -t pass_filtered_validated_variants.vcf.gz
                /opt/bcftools/bin/bcftools annotate -Oz -o validated_annotated_pipeline_variants.vcf.gz -a pass_filtered_validated_variants.vcf.gz -m 'VALIDATED' $PIPELINE_VCF
                /opt/bcftools/bin/bcftools index -t validated_annotated_pipeline_variants.vcf.gz
            elif [ "$#" -eq 1 ]; then
                cp $PIPELINE_VCF validated_annotated_pipeline_variants.vcf.gz
                cp $PIPELINE_VCF.tbi validated_annotated_pipeline_variants.vcf.gz.tbi
            else
                exit 1
            fi

baseCommand: ["/bin/bash", "annotate.sh"]

inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 1
        doc: "Each variant in this file that is also in the validated vcf file (if supplied) will be marked with a VALIDATED flag in its INFO field"
    validated_variants:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
        doc: "A vcf of previously discovered variants to be marked in the pipeline vcf; if not provided, this tool does nothing but rename the input vcf"

outputs:
    validated_annotated_vcf:
        type: File
        outputBinding:
            glob: "validated_annotated_pipeline_variants.vcf.gz"
        secondaryFiles: [.tbi]
