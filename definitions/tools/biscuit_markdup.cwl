#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit dedup"
baseCommand: ["/usr/bin/biscuit", "markdup"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite"

arguments:
    ["/dev/stdout",
    { shellQuote: false, valueFrom: "|" },
    "/usr/bin/sambamba", "sort", "-t", $(runtime.cores), "-m", "8G", "-o", "$(runtime.outdir)/markdup.bam", "/dev/stdin"
    ]
inputs:
    bam:
        type: File
        inputBinding:
            position: -1
outputs:
    markdup_bam:
        type: File
        outputBinding:
            glob: "markdup.bam"
