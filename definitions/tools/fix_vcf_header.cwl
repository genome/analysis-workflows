#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "fix vcf header"

baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "FixVcfHeader"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/$(inputs.vcf.basename).gz }]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "broadinstitute/picard:2.23.6"
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "INPUT="
outputs:
    fixed_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.vcf.basename).gz"
        secondaryFiles: [.tbi]
