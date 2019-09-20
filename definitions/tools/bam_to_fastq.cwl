#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Picard: BAM to FASTQ"
baseCommand: ["/bin/bash","bam_to_fastq.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'bam_to_fastq.sh'
        entry: |
            set -eou pipefail

            bam=$1
            paired=$2            
            if [[ "$paired" == "true" ]];then
                /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I=$bam F=read1.fastq F2=read2.fastq
            else
                /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I=$bam F=read1.fastq
            fi

inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    paired_end: 
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
        inputBinding:
            position: 2
outputs:
    fastqs:
        type: File[]
        outputBinding:
            glob: "read*.fastq"
