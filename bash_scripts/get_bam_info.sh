#!/bin/bash

#Takes 1 positional argument: gms analysis project ID <analysis_ID>.
#Finds model_ids based of anlaysis project number. Creates gvcfs_hg38.list and cromwell.options
#for a finsihed model.

analysis_ID=$1
outDIR=/gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"
outfile=/gscmnt/gc2698/jin810/jointcalling/"$analysis_ID"/bam_info_.txt
model_ids=$(genome instrument-data --show=id --filter "sample.name:$analysis_ID")
DIR=/gscmnt/gc2698/jin810/model_data

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

