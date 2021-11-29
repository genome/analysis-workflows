#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "sequence (bam or fastqs) to trimmed fastqs"
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
    unzip_fastqs:
        type: boolean?
    
outputs:
    fastqs:
        type: File[]
        outputSource: trim_fastq/fastqs
    fastq1:
         type: File
         outputSource: trim_fastq/fastq1
    fastq2:
         type: File
         outputSource: trim_fastq/fastq2

steps:
    sequence_to_fastq:
        run: ../tools/sequence_to_fastq.cwl
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
            unzip_fastqs: unzip_fastqs
        out:
            [fastq1, fastq2]
    trim_fastq:
        run: ../tools/trim_fastq.cwl
        in:
            reads1: sequence_to_fastq/fastq1
            reads2: sequence_to_fastq/fastq2
            adapters: adapters
            adapter_trim_end: adapter_trim_end
            adapter_min_overlap: adapter_min_overlap
            max_uncalled: max_uncalled
            min_readlength: min_readlength
        out:
            [fastqs, fastq1, fastq2]
    
