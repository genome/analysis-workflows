#/bin/bash
#Runs exome GATK4 joint-calling and genotyping using hg38 and merged .g.vcf files. 
#Takes 1 input. Analysis ID as postiional argument 1.

#Input
analysis_ID=$1
DIR=/gscmnt/gc2698/jin810/jointcalling/$analysis_ID/merged_gvcfs
sample_array=$()
interval_list=/gscmnt/gc2698/jin810/references/jin_lab_exome_targets.interval_list

#/gscmnt/gc2698/jin810/bash_scripts/get_AnP_info.sh $analysis_ID &&

#Place each sample name into array without newline return (this was causing problems using only cut -f1 into an array variable)
for i in $( cut -f1 "$DIR"/gvcfs.list )
do
    echo $i
    sample_array+=($i)
    #echo file is ${file_array[$count]}
    #let "count+=1"
done

#Block iterates through each sample and associated merged.vcf.gz file to Select Variants from the given interval list
count=1
for n in $( cut -f2 "$DIR"/gvcfs.list )
do
    echo count is "$count"
    sample=${sample_array[$count]}
    echo sample is "$sample" and file is "$n"
    bsub -oo "$DIR/$sample"_output.txt -eo "$DIR/$sample"_error.txt -q research-hpc -a 'docker(broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2)' \
-R "select[mem>16000] rusage[mem=15000]" /usr/bin/java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.1.7.0-local.jar SelectVariants -V "$n" -L "$interval_list" -R /gscmnt/gc2560/core/model_data/ref_build_aligner_index_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/aligner-index-blade18-1-1.gsc.wustl.edu-tmooney-331-75fff591a14f4f7c910247fc39c4ea7f/bwamem/0_7_15/all_sequences.fa --remove-unused-alternates --preserve-alleles \
-O "$DIR/$sample"_intervals_subset.vcf.gz
    let "count+=1"  
done

