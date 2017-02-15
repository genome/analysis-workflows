#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "strelka 2.7.1"
baseCommand: "/usr/bin/docker_helper.pl"
requirements:
    DockerRequirement:
        dockerPull: "mgibio/strelka-cwl:2.7.1"
    ResourceRequirement:
        coresMin: 8
arguments:
    [ { valueFrom: $(runtime.cores), position: 1 },
      { valueFrom: $(runtime.outdir), position: 2 }]
inputs:
    tumor_bam:
        type: File
        inputBinding:
            prefix: '--tumorBam='
            separate: false
            position: 3
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        inputBinding:
            prefix: '--normalBam='
            separate: false
            position: 4
        secondaryFiles: [^.bai]
    reference:
        type: File
        inputBinding:
            prefix: '--referenceFasta='
            separate: false
            position: 5
        secondaryFiles: [.fai]
    exome_mode:
        type: boolean
        inputBinding:
            prefix: '--exome'
            position: 6
outputs:
     indels:
         type: File
         outputBinding:
             glob: "results/variants/somatic.indels.vcf.gz"
     snvs:
         type: File
         outputBinding:
             glob: "results/variants/somatic.snvs.vcf.gz"

