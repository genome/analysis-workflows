#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "bam to trimmed fastqs"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml

inputs:
    unaligned:
        type: ../types/sequence_data.yml#sequence_data
    adapters:
        type: File
    adapter_trim_end:
        type: string
    adapter_min_overlap:
        type: int
    max_uncalled:
        type: int
    min_readlength:
        type: int
    
outputs:
    fastqs:
        type: File[]
        outputSource: trim_fastq/fastqs
    fastq_1:
         type: File
         outputSource: trim_fastq/fastq_1
    fastq_2:
         type: File
         outputSource: trim_fastq/fastq_2

steps:
    bam_to_fastq:
        run: ../tools/sequence_to_fastq_rna.cwl
        in: 
            bam:
                source: unaligned
                valueFrom: "$(self.sequence.hasOwnProperty('bam')? self.sequence.bam : null)"
            fastq1:
                source: unaligned
                valueFrom: "$(self.sequence.hasOwnProperty('fastq1')? self.sequence.fastq1 : null)"
            fastq2:
                source: unaligned
                valueFrom: "$(self.sequence.hasOwnProperty('fastq2')? self.sequence.fastq2 : null)"
        out:
            [fastqW1, fastqW2]
    trim_fastq:
        run: ../tools/trim_fastq.cwl
        in:
            reads1: bam_to_fastq/fastqW1
            reads2: bam_to_fastq/fastqW2
            adapters: adapters
            adapter_trim_end: adapter_trim_end
            adapter_min_overlap: adapter_min_overlap
            max_uncalled: max_uncalled
            min_readlength: min_readlength
        out:
            [fastqs, fastq_1, fastq_2]
    
