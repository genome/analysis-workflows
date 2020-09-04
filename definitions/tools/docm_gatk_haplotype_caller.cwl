#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HaplotypeCaller (GATK4)"
baseCommand: ["/bin/bash", "docm_haplotypeCaller.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.2.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'docm_haplotypeCaller.sh'
        entry: |
             set -o pipefail
             set -o errexit

             # Running haplotype caller using the newly created interval list
             if [[ "$#" == 5 ]];then # If normal_cram is passed.
                 # explicitly capturing variables
                 reference=$1
                 normal_cram=$2
                 tumor_cram=$3
                 docm_vcf=$4
                 interval_list=$5
                 # Chaning the interval_list to a new docm_interval_list that spans the docm regions by 200bp
                 cat $interval_list | grep '^@' > docm.interval_list # Extracting the header from the interval_list
                 zcat $docm_vcf | grep ^chr | awk '{FS = "\t";OFS = "\t";print $1,$2-100,$2+100,"+",$1"_"$2-100"_"$2+100}' >> docm.interval_list # Extracting the docm regions with a 100bp flanking region on both directions
                 /gatk/gatk HaplotypeCaller --java-options "-Xmx8g" -R $reference -I $normal_cram -I $tumor_cram --alleles $docm_vcf -L docm.interval_list --genotyping-mode GENOTYPE_GIVEN_ALLELES -O docm_raw_variants.vcf
             else # If normal_cram is not passed
                 reference=$1
                 tumor_cram=$2
                 docm_vcf=$3
                 interval_list=$4
                 # Chaning the interval_list to a new docm_interval_list that spans the docm regions by 200bp
                 cat $interval_list | grep '^@' > docm.interval_list # Extracting the header from the interval_list
                 zcat $docm_vcf | grep ^chr | awk '{FS = "\t";OFS = "\t";print $1,$2-100,$2+100,"+",$1"_"$2-100"_"$2+100}' >> docm.interval_list # Extracting the docm regions with a 100bp flanking region on both directions
                 /gatk/gatk HaplotypeCaller --java-options "-Xmx8g" -R $reference -I $tumor_cram --alleles $docm_vcf -L docm.interval_list --genotyping-mode GENOTYPE_GIVEN_ALLELES -O docm_raw_variants.vcf
             fi

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            position: 1
    normal_cram:
        type: File?
        inputBinding:
            position: 2
        secondaryFiles: [^.crai]
    cram:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [^.crai]
    docm_vcf:
        type: File
        inputBinding:
            position: 4
        secondaryFiles: [.tbi]
    interval_list:
        type: File
        inputBinding:
            position: 5
outputs:
    docm_raw_variants:
        type: File
        outputBinding:
            glob: "docm_raw_variants.vcf"
