#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Trim FASTQ (flexbar)"
baseCommand: ['/bin/bash','trim_fastq.sh']
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      tmpdirMin: 25000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite:v1.4"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'trim_fastq.sh'
        entry: |
            set -eou pipefail

            paired="false"
            fastq2=""

            adapters="$1"
            adapter_trim_end="$2"
            adapter_min_overlap="$3"
            max_uncalled="$4"
            min_readlength="$5"
            target="$6"
            threads="$7"
            fastq1="$8"
            if [[ $# -gt 8 ]];then
                fastq2=$9
                paired="true"
            fi
            if [[ $# -gt 9 ]];then
               echo "ERROR: more than two fastq files provided to flexbar"
               exit 1
            fi            


            if [[ "$paired" == "true" ]];then                
                /opt/flexbar/flexbar --adapters $adapters --adapter-trim-end $adapter_trim_end --adapter-min-overlap $adapter_min_overlap --max-uncalled $max_uncalled --min-read-length $min_readlength --target $target --threads $threads --reads $fastq1 --reads2 $fastq2
            else
                /opt/flexbar/flexbar --adapters $adapters --adapter-trim-end $adapter_trim_end --adapter-min-overlap $adapter_min_overlap --max-uncalled $max_uncalled --min-read-length $min_readlength --target $target --threads $threads --reads $fastq1
            fi
arguments: [
    {valueFrom: "$(runtime.outdir)/trimmed_read", position: 6},
    {valueFrom: "$(runtime.cores)", position: 7}
]

inputs:
    adapters:
        type: File
        inputBinding:
            position: 1
    adapter_trim_end:
        type: string
        inputBinding:
            position: 2
    adapter_min_overlap:
        type: int
        inputBinding:
            position: 3
    max_uncalled:
        type: int
        inputBinding:
            position: 4
    min_readlength:
        type: int
        inputBinding:
            position: 5
    fastqs:
        type: File[]
        inputBinding:
            position: 8
outputs:
    trimmed_fastqs:
        type: File[]
        outputBinding:
            glob: "trimmed_read*fastq"
