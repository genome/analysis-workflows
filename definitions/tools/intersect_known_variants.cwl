#!/usr/bin/env cwl-runner
 
cwlVersion: v1.0
class: CommandLineTool
label: "Intersect passing validated variants and passing pipeline variants for use in pvacseq"

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
                VALIDATED_VCF="$2"
                #filter the validated vcf to ensure there are only passing variants, then re-index
                /opt/bcftools/bin/bcftools view -f PASS -Oz -o pass_filtered_validated_variants.vcf.gz $VALIDATED_VCF
                /opt/bcftools/bin/bcftools index -t pass_filtered_validated_variants.vcf.gz
                #intersect the two vcfs; output will contain only passing variants
                #-n specifies that the output should contain only variants found in both files
                #-w results in a single output vcf containing the intersection
                #-p specifies the directory that will contain output files (vcf, index, and summary files)
                #-Oz specifies the output format as compressed
                /opt/bcftools/bin/bcftools isec -f PASS -n=2 -w1 -p validated -Oz $PIPELINE_VCF pass_filtered_validated_variants.vcf.gz
            elif [ "$#" -eq 1 ]; then
                mkdir validated 
                cp $PIPELINE_VCF validated/0000.vcf.gz
                cp $PIPELINE_VCF.tbi validated/0000.vcf.gz.tbi
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
        doc: "Pipeline variants to be intersected with validated variants, if the vcf is present"
    validated_variants:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
        doc: "A vcf of previously discovered variants; if not provided, this tool does nothing but rename the input vcf"

outputs:
    validated_and_pipeline_vcf:
        type: File
        outputBinding:
            glob: "validated/0000.vcf.gz"
        secondaryFiles: [.tbi]
