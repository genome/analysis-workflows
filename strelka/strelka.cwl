#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "strelka 1.0.15"
baseCommand: "/usr/bin/cwl_helper.sh"
requirements:
    DockerRequirement:
        dockerPull: "mgibio/strelka-cwl:1.0.15"
    ResourceRequirement:
        coresMin: 8
arguments:
    - { valueFrom: $(runtime.cores), position: 5 }
inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: .bai
    normal_bam:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: .bai
    reference:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: .fai
    config:
        type: File
        inputBinding:
            position: 4
outputs:
     all_indels:
         type: File
         outputBinding:
             glob: "output/results/all.somatic.indels.vcf"
     all_snvs:
         type: File
         outputBinding:
             glob: "output/results/all.somatic.snvs.vcf"
     passed_indels:
         type: File
         outputBinding:
             glob: "output/results/passed.somatic.indels.vcf"
     passed_snvs:
         type: File
         outputBinding:
             glob: "output/results/passed.somatic.snvs.vcf"


