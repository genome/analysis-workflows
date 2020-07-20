#!/bin/bash

#Takes 2 positional arguments: <input> A file containg a list of all .g.vcf chromosome files created by gatk HaploTypeCaller \
#<DIR> an directory specifiying where the output will be created.
#analysis_ID=$1
#DIR=$2
#samples=$( (ls "$DIR"/"analysis_ID"/*)

input=$1
DIR=$2
#Mounts local volume within docker image. This is necessary to be able to use data from the cluster!
export LSF_DOCKER_VOLUMES="$DIR:$DIR"
export LSF_DOCKER_NETWORK=host
export LSF_DOCKER_IPC=host


#Script for Picard MergeVcfs.This is a helper script for get_chrom_vcfs.
#java -jar /opt/picard/picard.jar MergeVcfs I="$DIR"/"analysis_ID"/"$input" O="$DIR"/"$analysis_ID"/merged_chroms.vcf.gz

/usr/bin/java -jar /opt/picard/picard.jar MergeVcfs I="$input" O="$DIR"/merged_chroms.vcf.gz
