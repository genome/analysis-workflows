#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit dedup"
baseCommand: ["/usr/bin/biscuit", "markdup"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 24000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite:v1.3"

arguments:
    ["/dev/stdout",
    { shellQuote: false, valueFrom: "|" },
    "/usr/bin/sambamba", "sort", "-t", $(runtime.cores), "-m", "15G", "-o", "$(runtime.outdir)/markdup.bam", "/dev/stdin"
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
