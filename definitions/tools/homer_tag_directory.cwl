#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Creating the homer tag directory'
baseCommand: ["/bin/bash", "homer_tag_directory.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "chrisamiller/homer"
    - class: ResourceRequirement
      ramMin: 32000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'homer_tag_directory.sh'
        entry: |
		set -o pipefail
		set -o errexit

		name=$1
		bam=$2
		reference_fasta=$3  #assumes that there's a .genome file next to the reference .fa
		#run in chrisamiller/homer container
		export PATH=$PATH:/opt/homer/bin
		echo "name:     $name"
		echo "bam:      $bam"
		echo "reffasta: $reference_fasta"
		if [[ ! -d $name/homer ]];then
  			mkdir -p $name/homer
		fi
		echo "creating tagDir"
		makeTagDirectory $name/homer $bam
		echo "tagdir2Bed"
		/gscuser/dspencer/programs/homer/bin/tagDir2bed.pl $name/homer/ | awk '$2>0 && $3>$2' | awk -v n=$name '{print >n"/"n".bed.part."$1}'
		#low mem way to ensure contigs are ordered appropriately
		cut -f 1 $(dirname $reference_fasta)/$(basename $reference_fasta .fa).genome | while read i;do if [[ -e $name/$name.bed.part.$i ]];then sort -nk 2 $name/$name.bed.part.$i;fi;done | bgzip -c >$name/$name.bed.gz
		rm -f $name/$name.bed.part.*
		echo "tabix"
		tabix -p bed $name/$name.bed.gz            

inputs:
    tag_directory_name:
        type: string
        inputBinding:
            position: 1
    sam:
        type: File
        inputBinding:
            position: 2
    reference:
        type: string
        inputBinding:
            position: 3
outputs:
    tag_directory:
        type: Directory
        outputBinding:
            glob: "$(inputs.tag_directory_name)"
