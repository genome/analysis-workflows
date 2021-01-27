#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "adapter for sequence_to_biscuit_alignments"
doc: "Some workflow engines won't stage files in our nested structure, so parse it out here"
requirements:
    - class: InlineJavascriptRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
    - class: StepInputExpressionRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data
        doc: "the unaligned sequence data with readgroup information"
    trimming_options:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    reference_index:
        type: string
outputs:
    aligned_bam:
        type: File
        outputSource: biscuit_align/aligned_bam
steps:
     biscuit_align:
        run: ../tools/biscuit_align.cwl
        in:
            bam:
                source: sequence
                valueFrom: "$(self.sequence.hasOwnProperty('bam')? self.sequence.bam : null)"
            fastq1:
                source: sequence
                valueFrom: "$(self.sequence.hasOwnProperty('fastq1')? self.sequence.fastq1 : null)"
            fastq2:
                source: sequence
                valueFrom: "$(self.sequence.hasOwnProperty('fastq2')? self.sequence.fastq2 : null)"
            read_group:
                source: sequence
                valueFrom: $(self.readgroup)
            trimming_options: trimming_options
            reference_index: reference_index
        out:
            [aligned_bam]
