#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Sambamba: merge"
baseCommand: ["/bin/bash", "merge.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/bam-merge:0.1"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'merge.sh'
        entry: |
            #!/bin/bash
            set -o pipefail
            set -o errexit
            set -o nounset

            SORTED=false

            while getopts "t:m:n:s:" opt; do
                case "$opt" in
                    t)
                        NTHREADS="$OPTARG"
                        ;;
                    m)
                        MEM="$OPTARG"
                        ;;
                    n)
                        OUTFILENAME="$OPTARG"
                        ;;
                    s)
                        SORTED=true
                        ;;
                esac
            done

            BAMS=("${@:$OPTIND}")
            NUM_BAMS=${#BAMS[@]}

            #if there is only one bam, just copy it and index it
            if [[ $NUM_BAMS -eq 1 ]]; then
                cp "$BAMS" "$OUTFILENAME";
            else
                if [[ $SORTED == "true" ]];then
                    /usr/bin/sambamba merge -t "$NTHREADS" "$OUTFILENAME" "$BAMS"
                else #unsorted bams, use picard
                    cmd="java -jar -Xmx"
                    cmd+=$MEM
                    cmd+="m /opt/picard/picard.jar MergeSamFiles OUTPUT=$OUTFILENAME ASSUME_SORTED=true USE_THREADING=true SORT_ORDER=unsorted VALIDATION_STRINGENCY=LENIENT"
                    for i in "${BAMS[@]}";do 
                      cmd+=" INPUT=$i"
                    done
                    `$cmd`;
                fi
            fi
            if [[ $SORTED == "true" ]];then
              /usr/bin/sambamba index "$OUTFILENAME"
            fi

arguments: [
    "-t", "$(runtime.cores)",
    "-m", "$(runtime.ram)",
]
inputs:
    bams:
        type: File[]
        inputBinding:
            position: 3
    sorted:
        type: boolean?
        default: "false"
        inputBinding:
            prefix: "-s"
            position: 2
    name:
        type: string?
        default: "merged.bam"
        inputBinding:
            prefix: "-n"
            position: 1
outputs:
    merged_bam:
        type: File
        outputBinding:
            glob: "$(inputs.name)"
