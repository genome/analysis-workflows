#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "strelka 2.9.9"
baseCommand: ["/usr/bin/perl", "/usr/bin/docker_helper.pl"]
requirements:
    ResourceRequirement:
        coresMin: 4
        ramMin: 4000
    DockerRequirement:
      dockerPull: "mgibio/strelka-cwl:2.9.9"
arguments:
    [ { valueFrom: $(inputs.cpu_reserved), position: 1 },
      { valueFrom: $(runtime.outdir), position: 2 }]
inputs:
    tumor_cram:
        type: File
        inputBinding:
            prefix: '--tumorBam='
            separate: false
            position: 3
        secondaryFiles: [.crai,^.crai]
    normal_cram:
        type: File
        inputBinding:
            prefix: '--normalBam='
            separate: false
            position: 4
        secondaryFiles: [.crai,^.crai]
    reference:
        type: string
        inputBinding:
            prefix: '--referenceFasta='
            separate: false
            position: 5
    exome_mode:
        type: boolean
        inputBinding:
            prefix: '--exome'
            position: 6
    cpu_reserved:
        type: int?
outputs:
     indels:
         type: File
         outputBinding:
             glob: "results/variants/somatic.indels.vcf.gz"
     snvs:
         type: File
         outputBinding:
             glob: "results/variants/somatic.snvs.vcf.gz"

