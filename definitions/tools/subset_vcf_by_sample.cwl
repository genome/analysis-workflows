#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "subset a multi-sample vcf so that only a single sample is retained"
baseCommand: ["/bin/bash", "subset_vcf.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: 'zlskidmore/samtools:1.10'
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'subset_vcf.sh'
        entry: |
            set -eou pipefail

            vcf="$1" #csv list of vcfs
            sample="$2" #csv list of tumor sample names

            #die if the specified sample name doesn't exist
            sample_exists="FALSE"
            for i in `gunzip -c "$vcf" | grep "^#CHROM" | cut -f 10-`;do 
              if [[ "$i" == "$sample" ]];then
                sample_exists="TRUE"
              fi
            done
            if [[ $sample_exists == "FALSE" ]];then
              echo "ERROR: Sample $sample does not exist in the VCF provided";
              exit 1
            fi

            /usr/local/bin/bcftools view -s $sample $vcf | /usr/local/bin/bgzip -c >subset.vcf.gz;
            /usr/local/bin/tabix -p vcf subset.vcf.gz;


inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
        doc: "multi-sample vcf"
    sample_name:
        type: string
        inputBinding:
            position: 2
        doc: "sample name to keep"
outputs:
    subset_vcf:
        type: File
        outputBinding:
            glob: "subset.vcf.gz"
        secondaryFiles: [".tbi"]
