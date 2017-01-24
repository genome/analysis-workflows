#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'revert to fastq'
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "SamToFastq"]
arguments:
    ["FASTQ=", { valueFrom: $(runtime.outdir)/revert_to_fastq.rFastq1 },
    "SECOND_END_FASTQ=", { valueFrom: $(runtime.outdir)/revert_to_fastq.rFastq2 }]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
            position: 1
outputs:
    fastq:
        type: File
        outputBinding:
            glob: "revert_to_fastq.rFastq1"
    second_end_fastq:
        type: File
        outputBinding:
            glob: "revert_to_fastq.rFastq2"
