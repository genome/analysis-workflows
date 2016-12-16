#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "align with bwa_mem and tag"

baseCommand: ["/bin/bash", "helper.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
    - class: ResourceRequirement
      coresMin: 8
    - class: InitialWorkDirRequirement
      listing:
          - entryname: 'helper.sh'
            entry: |
                /usr/local/bin/bwa mem -K 100000000 -t $5 -Y -R "$1" $2 $3 $4 | /usr/local/bin/samblaster -a --addMateTags | /usr/local/bin/samtools view -b -S /dev/stdin
stdout: "refAlign.bam"
arguments:
    - position: 5
      valueFrom: $(runtime.cores)
inputs:
    readgroup:
        type: string
        inputBinding:
            position: 1
    reference:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac, .alt]
    fastq:
        type: File
        inputBinding:
            position: 3
    fastq2:
        type: File
        inputBinding:
            position: 4
outputs:
    aligned_bam:
            type: stdout
