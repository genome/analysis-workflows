#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: Quant"
baseCommand: ["/bin/bash","kallisto.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'kallisto.sh'
        entry: |
            set -eou pipefail

            ncores="$1"
            index="$2"
            strand="$3"
            paired="$4"

            #fragment params are ignored and inferred from the data for paired-end
            fragment_length="$5"
            fragment_length_stddev="$6"

            fastqs="$7"
            if [[ "$paired" == "true" ]];then
              /usr/bin/kallisto quant -t $ncores -b 100 --fusion -o kallisto -i $index $strand $fastqs
            else
              ##"Reasonable" values for fragment length and stddev are hardcoded here. If your input data is substantially different
              ##then you'll need to specify that
              /usr/bin/kallisto quant -t $ncores -b 100 --fusion -o kallisto -l $fragment_length -s $fragment_length_stddev -i $index $strand --single $fastqs
            fi

arguments: [
    $(runtime.cores)
]
inputs:
    kallisto_index:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    strand:
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.strand) {
                        if (inputs.strand == 'first') {
                            return ['--fr-stranded'];
                        } else if (inputs.strand == 'second') {
                            return ['--rf-stranded'];
                        } else {
                            return [];
                        }
                    } else {
                            return [];
                    }
                }
            position: 2
    paired_end:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
        inputBinding:
            position: 3
    fragment_length:
        type: int?
        default: 200
        inputBinding:
            position: 4
    fragment_length_stddev:
        type: int?
        default: 30
        inputBinding:
            position: 5
    fastqs:
        type:
            type: array
            items:
                type: array
                items: File
        inputBinding:
            position: 6
            itemSeparator: " "

outputs:
    expression_transcript_table:
        type: File
        outputBinding:
            glob: "kallisto/abundance.tsv"
    expression_transcript_h5:
        type: File
        outputBinding:
            glob: "kallisto/abundance.h5"
    fusion_evidence:
        type: File
        outputBinding:
            glob: "kallisto/fusion.txt"
