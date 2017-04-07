#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "strelka 2.7.1"
baseCommand: ["/usr/bin/perl", "/usr/bin/strelka_helper.pl"]
requirements:
    ResourceRequirement:
        coresMin: 8
        ramMin: 4000
arguments:
    [ { valueFrom: $(runtime.cores), position: 1 },
      { valueFrom: $(runtime.outdir), position: 2 }]
inputs:
    tumor_cram:
        type: File
        inputBinding:
            prefix: '--tumorBam='
            separate: false
            position: 3
        secondaryFiles: [.crai]
    normal_cram:
        type: File
        inputBinding:
            prefix: '--normalBam='
            separate: false
            position: 4
        secondaryFiles: [.crai]
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

