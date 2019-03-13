#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit pileup"
baseCommand: ["/bin/bash", "helper.sh"]
stdout: pileup.vcf.gz
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite:v1.3"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'helper.sh'
        entry: |
            set -eo pipefail

            cores=$1
            statsfile_output=$2
            reference_fasta=$3
            bam=$4
            
            /usr/bin/biscuit pileup -q $cores -w $statsfile_output $reference_fasta $bam | /opt/htslib/bin/bgzip
arguments: [
    { valueFrom: $(runtime.cores), position: -9 },
    { valueFrom: "pileup_stats.txt", position: -8 }
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
