#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HISAT2: align"
baseCommand: ["bash hisat.sh"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 16
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'hisat.sh'
        entry: |
            set -eou pipefail

            cores="$1"
            reference_index="$2"
            read_group_id="$3"
            read_group_fields"$4"
            paired="$5"
            strand="$6"
            #Shandle one or two fastq inputs
            fastq1="$7"
            if [[ $# -gt 7 ]];then
                fastq2=$8
            fi
            if [[ $# -gt 8 ]];then
               echo "ERROR: too many arguments - were more than two fastq files provided to HISAT?"
               exit 1
            fi

            #collapse readgroups into string
            rg=""
            for i in $read_group_fields;do 
                rg+="--rg $i "
            done

            strand_param=""

            if [[ $paired == "true" ]];then
                ##get the right strand parameter
                if [[ $strand == "first" ]];then
                    strand_param="--rna-strandness RF"
                elif [[ $strand == "second" ]];then
                    strand_param="--rna-strandness FR"
                fi
                /usr/bin/hisat2 -p "$cores" --dta -x reference --rg-id "$read_group_id" --rg "read_group_fields" -1 "$fastq1" -2 "$fastq2" "$strand_param" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o $(runtime.outdir)/aligned.bam /dev/stdin
            elif [[ "$paired" == "false" ]];then
                ##get the right strand parameter
                if [[ $strand == "first" ]];then
                    strand_param="--rna-strandness R"
                elif [[ $strand == "second" ]];then
                    strand_param="--rna-strandness F"
                fi  
                /usr/bin/hisat2 -p "$cores" --dta -x reference --rg-id "$read_group_id" --rg "read_group_fields" -1 "$fastq1" "$strand_param" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o $(runtime.outdir)/aligned.bam /dev/stdin
            fi
            
arguments: [
    {valueFrom: "$(runtime.cores)", position: 1}
]
inputs:
    reference_index:
        type: File
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
        inputBinding:
            position: 2
    fastqs:
        type: File[]
        inputBinding:
            position: -2
    read_group_id:
        type: string
        inputBinding:
            position: 3
    read_group_fields:  #how to collapse this?
        type:
            type: array
            items: string
            inputBinding:
        inputBinding:
            position: 4
    paired_end:
        type: boolean?
        default: true
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
    strand:
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
        inputBinding:
            position: 5
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
