#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "merge several somatic VCFs into a multi-sample VCF"
baseCommand: ["/bin/bash", "merge_somatic_vcfs.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: 'zlskidmore/samtools:1.10'
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'merge_somatic_vcfs.sh'
        entry: |
            set -eou pipefail

            vcfs="$1" #csv list of vcfs
            samples="$2" #csv list of tumor sample names
            normal_sample_name="$3"
            
            IFS=',' read -ra vcf_array <<< "$vcfs"
            IFS=',' read -ra sample_array <<< "$samples"
            
            if [ "${#vcf_array[@]}" != "${#sample_array[@]}" ]; then
                echo "ERROR: the number of vcfs must be equal to the number of sample names provided";
                exit 1
            fi
            
            # prep all the tumor vcfs for merging
            for i in ${!vcf_array[*]}; do
                samp=${sample_array[$i]};
                vcf=${vcf_array[$i]};
                #assumes the input VCF contains the samples in the order NORMAL TUMOR
                normal_name_in_vcf=`gunzip -c "$vcf" | grep "^#CHROM" | cut -f 10`;
                tumor_name_in_vcf=`gunzip -c "$vcf" | grep "^#CHROM" | cut -f 11`;
                #use index for name of temp files (instead of sample names) to avoid potential weird character issues
                /usr/local/bin/bcftools view -s "$tumor_name_in_vcf" "$vcf" | /usr/local/bin/bcftools reheader -s <(echo "$samp") | /usr/local/bin/bgzip -c >"$i.tumor_only.vcf.gz";
                /usr/local/bin/tabix -p vcf "$i.tumor_only.vcf.gz";
                
                #also create normal-only vcfs
                /usr/local/bin/bcftools view -s "$normal_name_in_vcf" "$vcf" | /usr/local/bin/bcftools reheader -s <(echo "$normal_sample_name") | /usr/local/bin/bgzip -c >"$i.normal_only.vcf.gz";
                /usr/local/bin/tabix -p vcf "$i.normal_only.vcf.gz";
            done
                
            #create a merged normal vcf that includes all variants from all samples - this isn't strictly necessary,
            #but gives us a more sane final vcf
            /usr/local/bin/bcftools concat -a -D *.normal_only.vcf.gz | /usr/local/bin/bgzip -c >merged.normal.vcf.gz
            /usr/local/bin/tabix -p vcf merged.normal.vcf.gz
            
            #finally, merge them all into one big VCF, keeping FILTER field as PASS if it is PASS in any sample
            cmd="/usr/local/bin/bcftools merge -F x merged.normal.vcf.gz"
            for i in ${!sample_array[*]};do
                cmd="$cmd $i.tumor_only.vcf.gz"
            done
            echo "Running: $cmd";
            `$cmd >merged.vcf`
            /usr/local/bin/bgzip merged.vcf
            /usr/local/bin/tabix -p vcf merged.vcf.gz

inputs:
    vcfs:
        type: File[]
        inputBinding:
            position: 1
            itemSeparator: ','
            separate: false
        doc: "input somatic bgzipped tabix indexed vcfs to merge. Expects sample order in vcf to be normal/tumor"
    sample_names:
        type: string[]
        inputBinding:
            position: 2
            itemSeparator: ','
            separate: false
        doc: "tumor sample names"
    normal_sample_name:
        type: string
        inputBinding:
            position: 3
        doc: "sample name of the shared normal sample"

outputs:
    merged_vcf:
        type: File
        outputBinding:
            glob: "merged.vcf.gz"
        secondaryFiles: [".tbi"]
