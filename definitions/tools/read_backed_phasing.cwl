#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Read-backed phasing"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "ReadBackedPhasing"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/gatk-cwl:3.6.0
arguments:
    ["--assumeIdenticalSamples",
     "-L", { valueFrom: $(inputs.vcf) },
     "-o", { valueFrom: $(runtime.outdir)/phased.vcf }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
    vcf:
        type: File
        inputBinding:
            prefix: "-V"
            position: 3
        secondaryFiles: [".tbi"]
outputs:
    phased_vcf:
        type: File
        outputBinding:
            glob: "phased.vcf"
