#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HISAT2: align RNAseq data"
baseCommand: ["/bin/bash","hisat.sh"]
requirements:
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

            while getopts "c:r:i:g:p:s:f:" opt; do
                case "$opt" in
                    c)
                        cores="$OPTARG"
                        ;;
                    r)
                        reference_index="$OPTARG"
                        ;;
                    i)
                        read_group_id="$OPTARG"
                        ;;
                    g)
                        read_group_fields="$OPTARG"
                        ;;
                    p)
                        paired="$OPTARG"
                        ;;
                    s)
                        strand="$OPTARG"
                        ;;
                    f)
                        fastq_string="$OPTARG"
                        ;;
                esac
            done


            #collapse readgroups into separastring
            rg=""
            for i in $read_group_fields;do
                rg+="--rg $i "
            done

            if [[ "$paired" == "true" ]];then
                ##get the right strand parameter
                if [[ "$strand" == "first" ]];then
                    strand_param="--rna-strandness RF"
                elif [[ "$strand" == "second" ]];then
                    strand_param="--rna-strandness FR"
                fi
                #split fastqs
                IFS=',' read -ra FASTQS <<< "$fastq_string";
                /usr/bin/hisat2 -p "$cores" --dta -x "$reference_index" --rg-id "$read_group_id" "$rg" -1 "${FASTQS[0]}" -2 "${FASTQS[1]}" "$strand_param" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o aligned.bam /dev/stdin

            elif [[ "$paired" == "false" ]];then
                ##get the right strand parameter
                if [[ "$strand" == "first" ]];then
                    strand_param="--rna-strandness R"
                elif [[ "$strand" == "second" ]];then
                    strand_param="--rna-strandness F"
                fi
                #assumes only one fastq passed in
                /usr/bin/hisat2 -p "$cores" --dta -x "$reference" --rg-id "$read_group_id" "$rg" -1 "$fastq_string" "$strand_param" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o aligned.bam /dev/stdin
            fi

arguments: [
    {valueFrom: "$(runtime.cores)", position: 1, prefix: "-c"}
]
inputs:
    reference_index:
        type: File
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
        inputBinding:
            position: 2
            prefix: "-r"
    fastqs:
        type: File[]
        inputBinding:
            position: 7
            prefix: "-f"
            itemSeparator: ","
    read_group_id:
        type: string
        inputBinding:
            position: 3
            prefix: "-i"
    read_group_fields:
        type:
            type: array
            items: string
        inputBinding:
            position: 4
            prefix: "-g"
            itemSeparator: " "
    paired_end:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
        inputBinding:
            position: 5
            prefix: "-p"
    strand:
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
        inputBinding:
            position: 6
            prefix: "-s"
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
