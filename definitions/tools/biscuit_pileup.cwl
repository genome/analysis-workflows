#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit pileup"
baseCommand: ["/bin/bash", "helper.sh"]
stdout: pileup.vcf.gz
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/biscuit:0.3.8"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'helper.sh'
        entry: |
            set -eo pipefail

            cores=$1
            reference_fasta="$2"
            bam="$3"
            
            /usr/bin/biscuit pileup -q $cores -w pileup_stats.txt  "$reference_fasta" "$bam" | /opt/htslib/bin/bgzip
arguments: [
    { valueFrom: $(runtime.cores), position: -9 },
]
inputs:
    bam:
        type: File
        inputBinding:
            position: -1
    reference: 
        type: string
        inputBinding:
            position: -2
outputs:
    vcf:
        type: stdout
