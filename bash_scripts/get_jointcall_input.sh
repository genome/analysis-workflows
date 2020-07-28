#!/bin/bash

#Takes 1 positional argument: gms analysis project ID <analysis_ID>.
#Finds model_ids based of anlaysis project number. Creates gvcfs_hg38.list and cromwell.options
#for a finsihed model.

analysis_ID=$1
outDIR=/gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"
outfile=/gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"/gvcfs_hg38.list
model_ids=$(genome model list --show=id --filter "analysis_project.id=$analysis_ID")
DIR=/gscmnt/gc2698/jin810/model_data
logs_dir="$outDIR"/logs
config_template=/gscmnt/gc2698/jin810/cromwell/cromwell_lsf_config_template.conf
config_file="$outDIR"/cromwell.conf
json=/gscmnt/gc2698/jin810/analysis-workflows/cromwell_wdl/jointgt_GATK4_exome_hg38_inputs.json
jt_runs="$outDIR"/jt_runs.txt
callsDIR="$outDIR"/cromwell-executions/JointGenotyping

#Main loop that iterates through model ids.

#Checks if directory already exists. Creates new if it does not.
if [ -d "$outDIR" ]
then
    echo "directory already made for analysis project"
else
    echo "Making directory for $outDIR"
    mkdir -p "$outDIR"
fi

#Makes new .list file for all merged gvcf files
if [ -f "$outfile" ]
then
    echo "gvcf.list already exists...deleting"
    rm "$outfile"
else
    echo "Making new gvcf.list"
fi
touch "$outfile"

#Copy cromwell java config file and modify for project specific outputs
cp "$config_template" "$outDIR"
mv "$outDIR/cromwell_lsf_config_template.conf" "$config_file"
sed -i "s+/gscmnt/gc2764/cad/jgarza/logs/cromwell-%J.err+$logs_dir/cromwell-%J.err+g" "$config_file" #error logs directory
sed -i "s+/gscmnt/gc2764/cad/jgarza/tmp/cromwell-executions+$outDIR/cromwell-executions+g" "$config_file" #root cromwell working directory
sed -i "s+/gscmnt/gc2764/cad/jgarza/logs/cromwell-workflow-logs+$logs_dir/cromwell-workflow-logs+g" "$config_file" #workflow logs directory

#List previous joint genotyping runs into a file
#WORK IN PROGRESS
#if [-d "$callsDIR"]
#then
    #if [-f "$jt_runs"]
    #then
        #echo "Previous joint genotyping run txt found!"
        #rm "$jt_runs"
    #else
        #echo "No previous jt_runs.txt found!"
#ls "$outDIR"/cromwell-executions > "$jt_runs"

# Copy cromwell json file and modify to incorporate output path
#WORK IN PROGRESS
#echo Current Joint Genotyping
#cp "$json" "$outDIR/"
#sed -i "4 i $outDIR" "$outDIR"/jointgt_GATK4_exome_hg38_inputs.json
#sed -i "4 i $jt_runs" "$outDIR"/jointgt_GATK4_exome_hg38_inputs.json

for i in $model_ids
do
    #Filters headers of genome model list by length of characters. ID's are about 32 characters in length.
    length=${#i}
    #echo id is "$i" and length is "$length"
    if [ $length -gt 5 ]
    then
        model_sub=( $(genome model list --show=subject --filter "genome_model_id=$i") ) #Saves output as array
        echo subject is ${model_sub[2]}

        #Block below finds builds for all models of an analysis project, gathers .gvcfs, and submits MergeVcf on the gvcfs for each sample (build)
        #Produces an intermediate file that contains a list of all the .gvcf paths named chrom_gvcfs.list.
        build=$(ls "$DIR/$i"/)
        for x in $build
        do
	    outFinal="$outDIR" #Final directory location. Output is now directed here.
            
            #Checks for *merged.vcf.gz in output directory of model_data/MODEL/BUILD/results
            if [ -f "$DIR/$i/$x/results/"*merged.vcf.gz ]
            then
                vcf=$(ls "$DIR/$i/$x/results/"*merged.vcf.gz)
                tbi=$(ls "$DIR/$i/$x/results/"*merged.vcf.gz.tbi)
                echo -e "${model_sub[2]}\t$vcf\t$tbi\t$DIR/$i/$x/results"  >> "$outfile"
            else
                echo "merged.vcf.gz file not found for "$x""
            fi

        done
    fi
done

#Checks for a logs directory. Creates one if not found.
if [ -d "$logs_dir" ]
then
    echo logs directory already created
    rm "$logs_dir"/*
else
    mkdir "$logs_dir"
    echo created logs directory
fi

#Create cromwell.options file
echo -e "{\n\t"'"final_workflow_outputs_dir"': '"'"$outDIR"'"', \
"\n\t"'"final_workflow_log_dir"': '"'"$logs_dir"'"', \
"\n\t"'"final_call_logs_dir"': '"'"$logs_dir"'"' \
"\n}" > "$outDIR"/cromwell.options
