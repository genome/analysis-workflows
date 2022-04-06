#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "runs how_are_we_stranded_here to determine RNAseq data strandedness"

doc: |
  Uses how_are_we_stranded_here, a python package for testing strandedness. Runs Kallisto and Rseqc (infer-experiment-py) to
  to check which direction reads align once mapped in transcripts. It first creates a Kallisto index (or uses a pre-made index)
  of your organism's transcriptome. It then maps a small subset of reads (default 200000) to the transcriptome and uses Kallisto's
  --genomebam argument to project pseudoalignments to the genome sorted BAM file. (Currently only Kallisto version 0.44.0 works well with
  how_are_we_stranded_here.) It finally runs RSeQC's infer_experiment.py to check which direction reads from the first and second pairs
  are aligned in relation to the transcript strand, and provides output with the likely strandedness of your data.

baseCommand: ['check_strandedness']
arguments: [
    "--print_commands"
]

stdout: "$(inputs.reads1.nameroot).strandness_check.txt"

requirements:
    - class: ResourceRequirement
      ramMin: 4000
      tmpdirMin: 5000
      coresMin: 1
    - class: DockerRequirement
      dockerPull: "mgibio/checkstrandedness:v1"
inputs:
    reference_annotation:
        type: File
        inputBinding:
            prefix: --gtf
    kallisto_index:
        type: File
        inputBinding:
            prefix: --kallisto_index
    cdna_fasta:
        type: File
        inputBinding:
            prefix: --transcripts
    reads1:
        type: File
        inputBinding:
            prefix: "--reads_1"
    reads2:
        type: File
        inputBinding:
            prefix: "--reads_2"
outputs:
    strandedness_check:
        type: stdout
