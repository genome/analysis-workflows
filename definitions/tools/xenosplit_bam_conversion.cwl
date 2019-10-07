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
            
            /opt/samtools/bin/samtools view -h -F 256 -F 2048 $1 > graftbam.bam
            /opt/samtools/bin/samtools sort -n -o graftbam_accepted.bam graftbam.bam
            /opt/samtools/bin/samtools index graftbam_accepted.bam

            /opt/samtools/bin/samtools view -h -F 256 -F 2048 $2 > hostbam.bam
            /opt/samtools/bin/samtools sort -n -o hostbam_accepted.bam hostbam.bam
            /opt/samtools/bin/samtools index hostbam_accepted.bam

inputs:
    graftbam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [.bai]
    hostbam:
        type: File?
        inputBinding:
            position: 2
        secondaryFiles: [.bai]

outputs:
    graftbam_accepted:
        type: File
        outputBinding:
            glob: "graftbam_accepted.bam"
        secondaryFiles: [.bai]
    hostbam_accepted:
        type: File
        outputBinding:
            glob: "hostbam_accepted.bam"
        secondaryFiles: [.bai]
