#!/usr/bin/env cwl-runner
 
cwlVersion: v1.0
class: CommandLineTool
label: "Adds an INFO tag (PREVIOUSLY_DISCOVERED) flagging variants in the target vcf present in a known-variants file"

requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"
    - class: ResourceRequirement
      ramMin: 8000
    - class: StepInputExpressionRequirement

baseCommand: ["/opt/bcftools/bin/bcftools", "annotate"]
arguments:
    [ "-Oz", "-o", "known_variants_filtered.vcf.gz" ]  

inputs:
    known_variants:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            position: 1
            valueFrom: |
                ${
                    return [ '-a', self.path, '-m', 'PREVIOUSLY_DISCOVERED' ];
                }
        doc: "A vcf of previously discovered variants to be marked in the second input vcf; if not provided, this tool does nothing but rename the second input vcf"
    vcf:
        type: File
        secondaryFiles: [.tbi]
        inputBinding:
            position: 2
        doc: "Each variant in this file that is also in the above file (if supplied) will be marked with a PREVIOUSLY_DISCOVERED flag in its INFO field"
outputs:
    known_filtered:
        type: File
        outputBinding:
            glob: "known_variants_filtered.vcf.gz"
