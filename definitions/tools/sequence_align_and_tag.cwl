#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "align with bwa_mem and tag"
doc: "Due to workflow runner limitations, use sequence_align_and_tag_adapter.cwl subworkflow to call this"
baseCommand: ["/bin/bash", "sequence_alignment_helper.sh"]
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 20000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/alignment_helper-cwl:1.1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'sequence_alignment_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit
            set -o nounset

            RUN_TRIMMING="false"

            while getopts "b:?1:?2:?g:r:n:t:?o:?" opt; do
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

            if [[ "$MODE" == 'fastq' ]]; then
                if [[ "$RUN_TRIMMING" == 'false' ]]; then
                    /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -R "$READGROUP" "$REFERENCE" "$FASTQ1" "$FASTQ2" | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
                else
                    /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads "$FASTQ1" --reads2 "$FASTQ2" --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
                      | /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
                fi
            fi
            if [[ "$MODE" == 'bam' ]]; then
                if [[ "$RUN_TRIMMING" == 'false' ]]; then
                    /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout | /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
                else
                   /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout \
                     | /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads - --interleaved --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
                     | /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
                fi
            fi
stdout: "refAlign.bam"
arguments:
    - valueFrom: $(runtime.cores)
      position: 5
      prefix: '-n'
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
    readgroup:
        type: string
        inputBinding:
            prefix: '-g'
    reference:
        type:
            - string
            - File
        secondaryFiles: [.amb, .ann, .bwt, .pac, .sa]
        inputBinding:
            position: 4
            prefix: '-r'
        doc: 'bwa-indexed reference file'
    trimming:
        type:
          - ../types/trimming_options.yml#trimming_options
          - "null"
        inputBinding:
            valueFrom: $( ['-t', self.adapters.path, '-o', self.min_overlap] )
outputs:
    aligned_bam:
        type: stdout
