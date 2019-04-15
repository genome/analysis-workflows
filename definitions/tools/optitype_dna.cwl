#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Run optitype on dna samples'
baseCommand: ["/bin/bash", "/usr/bin/optitype_script.sh"]
arguments:
    [ { valueFrom: $(runtime.outdir) },
    { valueFrom: $(runtime.outdir) } ]
requirements:
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/immuno_tools-cwl:1.0.0"
inputs:
    optitype_name:
        type: string?
        default: "optitype"
        doc: "A prefix that will be used to name all files produced by the script"
        inputBinding:
            position: 1
    cram:
        type: File
        doc: "File to be HLA-typed"
        inputBinding:
            position: 2
outputs:
    optitype_tsv:
        type: File
        outputBinding:
            glob: $(inputs.optitype_name)_result.tsv
    optitype_plot:
        type: File
        outputBinding:
            glob: $(inputs.optitype_name)_coverage_plot.pdf
