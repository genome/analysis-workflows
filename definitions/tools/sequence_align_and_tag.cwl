#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "align with bwa_mem and tag"

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
    - class: InlineJavascriptRequirement
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
                /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" "$FASTQ1" "$FASTQ2" | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
            fi
            if [[ "$MODE" == 'bam' ]]; then
                /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I="$BAM" INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout | /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" /dev/stdin | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
            fi
stdout: "refAlign.bam"
arguments:
    - valueFrom: $(runtime.cores)
      position: 5
      prefix: '-n'
    - valueFrom: "$(inputs.unaligned.sequence.hasOwnProperty('bam')? inputs.unaligned.sequence.bam : null)"
      prefix: '-b'
    - valueFrom: "$(inputs.unaligned.sequence.hasOwnProperty('fastq1')? inputs.unaligned.sequence.fastq1 : null)"
      prefix: '-1'
    - valueFrom: "$(inputs.unaligned.sequence.hasOwnProperty('fastq2')? inputs.unaligned.sequence.fastq2 : null)"
      prefix: '-2'
    - valueFrom: $(inputs.unaligned.readgroup)
      prefix: '-g'
inputs:
    unaligned:
        type: ../types/sequence_data.yml#sequence_data
        doc: "the unaligned sequence data with readgroup information"
    reference:
        type: string
        inputBinding:
            position: 4
            prefix: '-r'
        doc: 'bwa-indexed reference file'
outputs:
    aligned_bam:
        type: stdout
