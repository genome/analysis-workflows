#!/usr/bin/env cwl-runner
 
cwlVersion: v1.0
class: CommandLineTool
label: "Adds an INFO tag (CLE) to variants in the target file present in a CLE-validated vcf file"

requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: ResourceRequirement
      ramMin: 8000

baseCommand: ["/opt/bcftools/bin/bcftools", "annotate"]
arguments:
    [ "-Oz", "-o", "cle_variants_flagged.vcf.gz" ]  

inputs:
    cle_variants:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            position: 1
            valueFrom: |
                ${
                    return [ '-a', inputs.cle_variants.path, '-m', 'CLE' ];
                }
    vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
outputs:
    cle_flagged:
        type: File
        outputBinding:
            glob: "cle_variants_flagged.vcf.gz"
