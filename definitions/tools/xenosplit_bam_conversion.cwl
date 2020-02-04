#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Converting input BAM files to Xenosplit friendly BAM files"

baseCommand: ["/bin/bash", "Xenosplit.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 20000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'Xenosplit.sh'
        entry: |
            set -o pipefail
            set -o errexit

            /opt/samtools/bin/samtools sort -o graftsorted.bam $1
            /opt/samtools/bin/samtools index graftsorted.bam
            /opt/samtools/bin/samtools view -h -F 256 -F 2048 graftsorted.bam > graftbam.bam
            /opt/samtools/bin/samtools sort -n -o graftbam_accepted.bam graftbam.bam

            /opt/samtools/bin/samtools sort -o hostsorted.bam $2
            /opt/samtools/bin/samtools index hostsorted.bam
            /opt/samtools/bin/samtools view -h -F 256 -F 2048 hostsorted.bam > hostbam.bam
            /opt/samtools/bin/samtools sort -n -o hostbam_accepted.bam hostbam.bam

inputs:
    graftbam:
        type: File
        inputBinding:
            position: 1
    hostbam:
        type: File?
        inputBinding:
            position: 2

outputs:
    graftbam_accepted:
        type: File
        outputBinding:
            glob: "graftbam_accepted.bam"
    hostbam_accepted:
        type: File
        outputBinding:
            glob: "hostbam_accepted.bam"
