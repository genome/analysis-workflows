#!/usr/bin/env cwl-runner
 
cwlVersion: v1.0
class: CommandLineTool
label: "Intersect passing cle variants and passing pipeline variants for use in pvacseq"

requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: ResourceRequirement
      ramMin: 8000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'intersect.sh'
        entry: |
            set -eou pipefail

            PIPELINE_VCF="$1"

            if [ "$#" -eq 2 ]; then
                CLE_VCF="$2"
                /opt/bcftools/bin/bcftools view -f PASS -Oz -o pass_filtered_cle_variants.vcf.gz $CLE_VCF
                /opt/bcftools/bin/bcftools index -t pass_filtered_cle_variants.vcf.gz
                /opt/bcftools/bin/bcftools isec -f PASS -n=2 -w1 -p cle -Oz $PIPELINE_VCF pass_filtered_cle_variants.vcf.gz
            elif [ "$#" -eq 1 ]; then
                mkdir cle 
                cp $PIPELINE_VCF cle/0000.vcf.gz
                cp $PIPELINE_VCF.tbi cle/0000.vcf.gz.tbi
            else
                exit 1
            fi

baseCommand: ["/bin/bash", "intersect.sh"]

inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 1
        doc: "Pipeline variants to be intersected with cle variants, if the vcf is present"
    cle_variants:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
        doc: "A vcf of previously discovered variants; if not provided, this tool does nothing but rename the input vcf"

outputs:
    cle_and_pipeline_vcf:
        type: File
        outputBinding:
            glob: "cle/0000.vcf.gz"
        secondaryFiles: [.tbi]
