#!/bin/bash

#Takes 2 positional argument: gms analysis project ID <analysis_ID> and <outDIR> for location of output dirs/files.
#Finds model_ids based of anlaysis project number. Outputs a bsub job to MergeVCFs for all chromosome files ( .g.vcf files)
#for a finsihed model.

analysis_ID=$1
outDIR=$2
model_ids=$(genome model list --show=id --filter "analysis_project.id=$analysis_ID")
DIR=/gscmnt/gc2698/jin810/model_data/

#Main loop that iterates through model ids.
#file_array=()
#path_array=()

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
	    outFinal="$outDIR"/"$analysis_ID"/${model_sub[2]}/"$i"/"$x"/ #Final directory location. Output is now directed here.
            data_path="$DIR/$i/$x/results"
            echo "$data_path"

            #Checks if directory already exists. Creates new if it does not.
            if [ -d "$outFinal"/${model_sub[2]}/"$i"/"$x"/ ]
            then
                echo "directory already made for subject"
            else
                echo "Making directory for ${model_sub[2]}"
                mkdir -p "$outFinal"/${model_sub[2]}/"$i"/"$x"/
            fi

            ls "$data_path"/*g.vcf.gz > "$outFinal"/chrom_gvcfs.list

            #file_array+=( "$outFinal"/chrom_gvcfs.list )
            #path_array+=( "$outFinal" )

            #Block submits a picard MergeVCFs job for each subject in analysis project.
            bsub -q research-hpc -R "select[mem>40000]rusage[mem=40000]" -a 'docker(mgibio/picard-cwl:2.18.1)' -oo "$outFinal"/mergevcfs_output.txt \
            -eo "$outFinal"/mergevcfs_error.txt /bin/bash ~/bash_scripts/merge_chrom_gvcfs.sh "$outFinal"/chrom_gvcfs.list "$outFinal"
        done
    fi
done
