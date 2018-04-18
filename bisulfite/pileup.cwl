#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit pileup"
baseCommand: ["/usr/bin/biscuit", "pileup"]
stdout: pileup.vcf.gz
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 4
arguments: [
    { valueFrom: "-q", position: -10 },
    { valueFrom: $(runtime.cores), position: -9 },
    { valueFrom: "-w", position: -8 },
    { valueFrom: "pileup_stats.txt", position: -7 },
    { shellQuote: false, valueFrom: "|" },
    "/bin/bgzip"
]
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "-i"
            position: -2
    reference: 
        type: string
        inputBinding:
            prefix: "-r"
            position: -1
outputs:
    vcf:
        type: stdout
