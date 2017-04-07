#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 somatic"
baseCommand: "/usr/bin/varscan_helper.sh"
requirements:
    - class: ResourceRequirement
      ramMin: 12000
inputs:
    tumor_cram:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.crai]
    normal_cram:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [^.crai]
    reference:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [.fai]
    roi_bed:
        type: File?
        inputBinding:
            position: 4
outputs:
    snvs:
        type: File
        outputBinding:
            glob: "output.snp.vcf"
    indels:
        type: File
        outputBinding:
            glob: "output.indel.vcf"
