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
          - $import: ../types/trimming_options.yml
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
    trimming:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    ref_fai:
        type: File
    ref_dict:
        type: File
    ref_amb:
        type: File
    ref_ann:
        type: File
    ref_bwt:
        type: File
    ref_pac:
        type: File
    ref_sa:
        type: File

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
            trimming: trimming
            ref_fai: ref_fai
            ref_dict: ref_dict
            ref_amb: ref_amb
            ref_ann: ref_ann
            ref_bwt: ref_bwt
            ref_pac: ref_pac
            ref_sa: ref_sa

        out:
            [aligned_bam]
