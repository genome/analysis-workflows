#!/bin/bash
#Runs exome GATK4 joint-calling and genotyping using hg38 and merged .g.vcf files. 
#Takes 1 input. Analysis ID as postiional argument 1.

#Input
analysis_ID=$1
main=/gscmnt/gc2698/jin810/cromwell/lek_joint_genotyping_jinlab_mod.wdl
json=/gscmnt/gc2698/jin810/cromwell/jointgt_GATK4_exome_hg38_inputs.json
options=/gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"/cromwell.options

#Main Script
cd /gscmnt/gc2698/jin810/jointcalling/"$analysis_ID" && \
bsub -q research-hpc -R "select[mem>5000] span[hosts=1] rusage[mem=3000]" -a "docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute0-24)" \
-oo /gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"/bsub_exec.out -eo /gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"/bsub_exec.err \
-J jc_"$analysis_ID" /usr/bin/java -Dsystem.input-read-limits.lines=18000000 -Dconfig.file="/gscmnt/gc2698/jin810/jointcalling/$analysis_ID/cromwell.conf" \
-jar /opt/cromwell.jar run "$main" \
-i "$json"
#-o "$options"
