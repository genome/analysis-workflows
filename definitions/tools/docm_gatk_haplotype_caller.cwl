#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HaplotypeCaller (GATK 3.6)"
baseCommand: ["/bin/bash", "docm_haplotypeCaller.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.6.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'docm_haplotypeCaller.sh'
        entry: |
             set -o pipefail
             set -o errexit

             # Extracting the regions from the interval list to match the DoCM VCF file
             cat $5 | grep '^@' > docm.interval_list
             zcat $4 | grep ^chr | awk '{FS = "\t";OFS = "\t";print $1,$2-100,$2+100,"+",$1"_"$2-100"_"$2+100}' >> docm.interval_list

             # Running haplotype caller using the newly created interval list
             if [[ "$2" eq "" ]];then
                 /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T HaplotypeCaller -R $1 -I $3 --alleles $4 -L docm.interval_list -gt_mode GENOTYPE_GIVEN_ALLELES -o docm_raw_variants.vcf # if normal_bam is not provided
             else
                 /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T HaplotypeCaller -R $1 -I $2 -I $3 --alleles $4 -L docm.interval_list -gt_mode GENOTYPE_GIVEN_ALLELES -o docm_raw_variants.vcf # if normal_bam is provided
             fi

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    normal_bam:
        type: File?
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [^.bai]
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 3
        secondaryFiles: [^.bai]
    docm_vcf:
        type: File
        inputBinding:
            prefix: "--alleles"
            position: 4
        secondaryFiles: [.tbi]
    interval_list:
        type: File
        inputBinding:
            prefix: "-L"
            position: 5
outputs:
    docm_raw_variants:
        type: File
        outputBinding:
            glob: "docm_raw_variants.vcf"
