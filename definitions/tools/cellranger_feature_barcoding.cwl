#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger Count with feature barcoding"

baseCommand: ["/bin/bash", "create_library_file.sh"]
arguments: [{ shellQuote: false, valueFrom: "&&" }, "/apps/cellranger-6.0.0/cellranger", "count", "--libraries=libraries.csv", "--id=$(inputs.run_name)", "--localcores=$(runtime.cores)", "--localmem=$(runtime.ram/1000)"]

requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/alex.paul/cellranger:6.0.0"
    - class: ResourceRequirement
      ramMin: 56000
      coresMin: 8
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "create_library_file.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail

          ANTIBODY_FASTQS="$1"
          ANTIBODY_SAMPLE_NAMES="$2"
          EXPRESSION_FASTQS="$3"
          EXPRESSION_SAMPLE_NAMES="$4"

          # add header
          echo "fastqs,sample,library_type" > libraries.csv

          # create arrays to iterate over
          IFS=',' read -ra ARRAY_ANTIBODY_FASTQS <<< "$ANTIBODY_FASTQS"
          IFS=',' read -ra ARRAY_ANTIBODY_NAMES <<< "$ANTIBODY_SAMPLE_NAMES"
          IFS=',' read -ra ARRAY_EXPRESSION_FASTQS <<< "$EXPRESSION_FASTQS"
          IFS=',' read -ra ARRAY_EXPRESSION_NAMES <<< "$EXPRESSION_SAMPLE_NAMES"

          if [ "${#ARRAY_ANTIBODY_FASTQS[@]}" != "${#ARRAY_ANTIBODY_NAMES[@]}" ] && [ "${#ARRAY_EXPRESSION_FASTQS[@]}" != "${#ARRAY_EXPRESSION_NAMES[@]}" ]; then
              echo "ERROR input sample names need to be equal to length of input fastq directories"
              exit 1
          fi

          # add antibody fastqs
          for i in "${!ARRAY_ANTIBODY_FASTQS[*]}"; do
              echo "${ARRAY_ANTIBODY_FASTQS[$i]},${ARRAY_ANTIBODY_NAMES[$i]},Antibody Capture" >> libraries.csv
          done

          # add expression fastqs
          for i in "${!ARRAY_EXPRESSION_FASTQS[*]}"; do
              echo "${ARRAY_EXPRESSION_FASTQS[$i]},${ARRAY_EXPRESSION_NAMES[$i]},Gene Expression" >> libraries.csv
          done
          exit 0

inputs:
    antibody_fastq_directory:
        type: Directory[]
        inputBinding:
           position: -4
           itemSeparator: ","
           separate: false
        doc: "Array of directories containing fastq files with antibody feature barcoding"
    antibody_sample_names:
        type: string[]
        inputBinding:
            position: -3
            itemSeparator: ","
            separate: false
        doc: "Sample name, must be same as name specified in sample sheet in previous mkfastq step"
    chemistry:
        type: string?
        inputBinding:
            prefix: --chemistry=
            position: 4
            separate: false
        default: "auto"
        doc: "Assay configuration used, default 'auto' should usually work without issue"
    expression_fastq_directory:
        type: Directory[]
        inputBinding:
           position: -2
           itemSeparator: ","
           separate: false
        doc: "Array of directories containing fastq files for gene expression"
    expression_sample_names:
        type: string[]
        inputBinding:
            position: -1
            itemSeparator: ","
            separate: false
        doc: "Sample name, must be same as name specified in sample sheet in previous mkfastq step"
    feature_reference_csv:
        type: File
        inputBinding:
            prefix: --feature-ref=
            position: 5
            separate: false
        doc: "Feature reference csv for the antibody barcodes"
    reference:
        type: Directory
        inputBinding:
            prefix: --transcriptome=
            position: 6
            separate: false
        doc: "Transcriptome reference compatible with input species and Cell Ranger"
    run_name:
        type: string
        doc: "Used to generate the run id"

outputs:
    libraries_csv:
        type: File
        outputBinding:
            glob: "libraries.csv"
    out_dir:
        type: Directory
        outputBinding:
            glob: "$(inputs.run_name)/outs/"

