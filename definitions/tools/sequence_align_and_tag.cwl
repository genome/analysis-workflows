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
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 20000
    - class: DockerRequirement
      dockerPull: "mgibio/alignment_helper-cwl:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'sequence_alignment_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit
            set -o nounset

            while getopts "b:?1:?2:?g:r:n:" opt; do
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
                esac
            done

            if [[ "$MODE" == 'fastq' ]]; then
                /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -R "$READGROUP" "$REFERENCE" "$FASTQ1" "$FASTQ2" | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
            fi
            if [[ "$MODE" == 'bam' ]]; then
                /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout | /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
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
outputs:
    aligned_bam:
        type: stdout
