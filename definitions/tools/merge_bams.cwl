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
                        SORTED="$OPTARG"
                        ;;
                esac
            done
            BAMS=${@:$OPTIND}
            NUM_BAMS=`echo "$# - $OPTIND + 1" | perl -nae 'print eval $_'` #can't use typical dollar/parens bash math

            #if there is only one bam, just copy it and index it
            if [[ $NUM_BAMS -eq 1 ]]; then
                cp "$BAMS" "$OUTFILENAME";
            else
                if [[ $SORTED == "true" ]];then
                    /usr/bin/sambamba merge -t "$NTHREADS" "$OUTFILENAME" "$BAMS"
                else #unsorted bams, use picard
                    cmd="java -jar -Xmx8g /opt/picard/picard.jar MergeSamFiles OUTPUT=$OUTFILENAME ASSUME_SORTED=true USE_THREADING=true SORT_ORDER=unsorted VALIDATION_STRINGENCY=LENIENT"
                    for i in $BAMS;do #this assumes no spaces in filenames, but a space in a filename is a space in one's soul
                      cmd="$cmd INPUT=$i"
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
    "-n", "$(inputs.name)",
    "-s", "$(inputs.sorted)"
]
inputs:
    bams:
        type: File[]
        inputBinding:
            position: 1
    sorted:
        type: string?
        default: "false"
    name:
        type: string?
        default: "output.bam"
outputs:
    merged_bam:
        type: File
        outputBinding:
            glob: "$(inputs.name)"
