#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Trim bisulfite reads and align usin biscuit"
baseCommand: ["/bin/bash","biscuit_trim_and_align.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 12
    - class: DockerRequirement
      dockerPull: "chrisamiller/biscuit:latest"
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'biscuit_trim_and_align.sh'
        entry: |
            set -o pipefail
            set -o errexit
            set -o nounset

            RUN_TRIMMING="false"

            while getopts "b:?1:?2:?g:r:n:d:t:?o:?" opt; do
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
                    g)
                        READGROUP="$OPTARG"
                        ;;
                    r)
                        REFERENCE="$OPTARG"
                        ;;
                    n)
                        NTHREADS="$OPTARG"
                        ;;
                    d)
                        OUTDIR="$OPTARG"
                        ;;
                    t)
                        RUN_TRIMMING="true"
                        TRIMMING_ADAPTERS="$OPTARG"
                        ;;
                    o)
                        RUN_TRIMMING="true"
                        TRIMMING_ADAPTER_MIN_OVERLAP="$OPTARG"
                        ;;
                esac
            done

            #reserve at least two cores for sorting
            #if there are too few for this, let the scheduler sort it out
            SORT_THREADS=2
            if [[ $NTHREADS -gt 3 ]];then
                NTHREADS=`expr $NTHREADS - $SORT_THREADS`
            fi
            #also reserve one core for trimming if it is specified
            if [[ "$RUN_TRIMMING" == "true" ]];then
                NTHREADS=`expr $NTHREADS - 1`
            fi

            if [[ "$MODE" == 'fastq' ]]; then
                if [[ "$RUN_TRIMMING" == 'false' ]]; then
                    /usr/bin/biscuit align -t "$NTHREADS" -M -R "$READGROUP" "$REFERENCE" "$FASTQ1" "$FASTQ2" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t "$SORT_THREADS" -m 8G -o "$OUTDIR/aligned.bam" /dev/stdin
                else
                    /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads "$FASTQ1" --reads2 "$FASTQ2" --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
                      | /usr/bin/biscuit align -t "$NTHREADS" -M -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t "$SORT_THREADS" -m 8G -o "$OUTDIR/aligned.bam" /dev/stdin
                fi
            fi
            if [[ "$MODE" == 'bam' ]]; then
                if [[ "$RUN_TRIMMING" == 'false' ]]; then
                    /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout | /usr/bin/biscuit align -t "$NTHREADS" -M -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t "$SORT_THREADS" -m 8G -o "$OUTDIR/aligned.bam" /dev/stdin
                else
                    /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout \
                    | /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads -r - \
                      | /usr/bin/biscuit align -t "$NTHREADS" -M -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t "$SORT_THREADS" -m 8G -o "$OUTDIR/aligned.bam" /dev/stdin
                fi
            fi

arguments: [
    { valueFrom: $(runtime.cores), position: -9, prefix: '-n' },
    { valueFrom: $(runtime.outdir), position: -8, prefix: '-d' },
]
inputs:
    reference_index:
        type: string
        inputBinding:
            position: -3
            prefix: '-r'
    bam:
        type: File?
        inputBinding:
            position: -5
            prefix: '-b'
    fastq1:
        type: File?
        inputBinding:
            position: -2
            prefix: '-1'
    fastq2:
        type: File?
        inputBinding:
            position: -1
            prefix: '-2'
    read_group:
        type: string
        inputBinding:
            position: -4
            prefix: '-g'
    trimming_options:
        type:
          - ../types/trimming_options.yml#trimming_options
          - "null"
        inputBinding:
            valueFrom: $( ['-t', self.adapters.path, '-o', self.min_overlap] )
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
