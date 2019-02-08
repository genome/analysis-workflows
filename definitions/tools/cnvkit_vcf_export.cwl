#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Convert default cnvkit .cns output to standard vcf format"

requirements:
    - class: DockerRequirement
      dockerPull: etal/cnvkit:0.9.5
    - class: ShellCommandRequirement
baseCommand: ["/usr/bin/python", "/usr/local/bin/cnvkit.py", "call"]
arguments: [
    { position: -1, valueFrom: $(inputs.male_reference), prefix: "-y" },
    "-o", "adjusted.tumor.cns",
    { shellQuote: false, valueFrom: "&&" },
    "/usr/bin/python", "/usr/local/bin/cnvkit.py", "export", "vcf", "adjusted.tumor.cns"
]
inputs:
    cns_file:
        type: File
        inputBinding:
            position: -2
    male_reference:
        type: boolean
        default: false
        inputBinding:
            position: 1
            prefix: "-y"
    cnr_file:
        type: File?
        inputBinding:
            position: 2
            prefix: "--cnr"
    output_name:
        type: string
        default: "cnvkit_output.vcf"
        inputBinding:
            position: 3
            prefix: "-o"
outputs:
    cnvkit_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_name)
