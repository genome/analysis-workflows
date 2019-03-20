#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Sort VCF"

baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/picard/picard.jar", "SortVcf"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/sorted.vcf }]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "I="
    reference_dict:
        type: File
        inputBinding:
            prefix: "SEQUENCE_DICTIONARY="
outputs:
    sorted_vcf:
        type: File
        outputBinding:
            glob: "sorted.vcf"
