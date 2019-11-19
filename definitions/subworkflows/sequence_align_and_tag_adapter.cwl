#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "adapter for sequence_align_and_tag"
doc: "Some workflow engines won't stage files in our nested structure, so parse it out here"
requirements:
    - class: InlineJavascriptRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
    - class: StepInputExpressionRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    unaligned:
        type: ../types/sequence_data.yml#sequence_data
        doc: "the unaligned sequence data with readgroup information"
    reference:
        type:
            - string
            - File
        secondaryFiles: [.amb, .ann, .bwt, .pac, .sa]
        doc: 'bwa-indexed reference file'
outputs:
    aligned_bam:
        type: File
        outputSource: align_and_tag/aligned_bam
steps:
    align_and_tag:
        run: ../tools/sequence_align_and_tag.cwl
        in:
            reference: reference
            bam:
                source: unaligned
                valueFrom: "$(self.sequence.hasOwnProperty('bam')? self.sequence.bam : null)"
            fastq1:
                source: unaligned
                valueFrom: "$(self.sequence.hasOwnProperty('fastq1')? self.sequence.fastq1 : null)"
            fastq2:
                source: unaligned
                valueFrom: "$(self.sequence.hasOwnProperty('fastq2')? self.sequence.fastq2 : null)"
            readgroup:
                source: unaligned
                valueFrom: $(self.readgroup)

        out:
            [aligned_bam]
