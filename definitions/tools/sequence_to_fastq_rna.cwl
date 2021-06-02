#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: "Picard: BAM to FASTQ"
baseCommand: ["/bin/bash","makefastqs.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 16000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq:1.0.0"

    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'makefastqs.sh'
        entry: |
            set -o pipefail
            set -o errexit
            set -o nounset
            while getopts ":b:?1:?2:d:?" opt; do
                case "$opt" in
                    b)
                        MODE=bam
                        BAM="$OPTARG"
                        ;;
                    1)
                        MODE=fastq
                        FASTQ1="$OPTARG"
                        ;;
                    2)  
                        MODE=fastq
                        FASTQ2="$OPTARG"
                        ;;
                    d)
                        OUTDIR="$OPTARG"
                        ;;
                esac
            done

            if [[ "$MODE" == 'fastq' ]]; then #must be fastq input
                # cp $FASTQ1 $OUTDIR/read.1.fastq
                # cp $FASTQ2 $OUTDIR/read.2.fastq
                gunzip -c $FASTQ1 > $OUTDIR/read1.fastq
                gunzip -c $FASTQ2 > $OUTDIR/read2.fastq
            else # then
                ##run samtofastq here, dumping to the same filenames
                ## input file is $BAM
                /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INCLUDE_NON_PF_READS=true F=$OUTDIR/read1.fastq F2=$OUTDIR/read2.fastq VALIDATION_STRINGENCY=SILENT
            fi
arguments: [
    {valueFrom: $(runtime.outdir), position: -6, prefix: '-d'}
]
inputs:
    bam:
        type: File?
        inputBinding:
            prefix: '-b'
    fastq1:
        type: File?
        inputBinding:
            prefix: '-1'
    fastq2:
        type: File?
        inputBinding:
            prefix: '-2'
outputs:
    fastqW1:
        type: File?
        outputBinding:
            glob: "read1.fastq"
    fastqW2:
        type: File?
        outputBinding:
            glob: "read2.fastq"
