#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard: BAM to FASTQ"
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/picard/picard.jar", "SamToFastq", "VALIDATION_STRINGENCY=SILENT"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq"
arguments: [ {valueFrom: "F=$(runtime.outdir)/read1.fastq"},
             {valueFrom: "F2=$(runtime.outdir)/read2.fastq"} ]
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
            separate: false
            position: 1
outputs:
    fastq1:
        type: File
        outputBinding:
            glob: "read1.fastq"
    fastq2:
        type: File
        outputBinding:
            glob: "read2.fastq"
