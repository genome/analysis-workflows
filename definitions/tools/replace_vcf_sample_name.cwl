#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/bin/bash", "rename_sample.sh"]
requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "rename_sample.sh"
        entry: |
          #!/bin/bash
          set -eou pipefail
          basen=`basename "$3"`
          basen="renamed.$basen"

          #escape spaces, otherwise bcftools will try to use them as a delimiter
          #triple backslash to escape within backticks and then again within sed
          old_name=`echo "$1" | sed 's/ /\\\ /g'`
          new_name=`echo "$2" | sed 's/ /\\\ /g'`

          echo "$old_name $new_name" > sample_update.txt
          /opt/bcftools/bin/bcftools reheader -s sample_update.txt -o "$basen" "$3"

inputs:
    input_vcf:
        type: File
        inputBinding:
            position: 3
        doc: "vcf file to filter"
    sample_to_replace:
        type: string
        inputBinding:
            position: 1
        doc: "Sample name to be replaced"
    new_sample_name:
        type: string
        inputBinding:
            position: 2
        doc: "Sample name to replace the other"

outputs:
    renamed_vcf:
        type: File
        outputBinding:
            glob: $("renamed." + inputs.input_vcf.basename)
