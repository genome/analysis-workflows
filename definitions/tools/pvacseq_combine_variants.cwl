#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Combine germline and somatic vcf for pVACseq phasing"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/gatk-cwl:3.6.0
arguments:
    ["--assumeIdenticalSamples",
     "-o", { valueFrom: $(runtime.outdir)/combined_somatic_plus_germline.vcf }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    germline_vcf:
        type: File
        inputBinding:
            prefix: "-V"
            position: 2
        secondaryFiles: [".tbi"]
    somatic_vcf:
        type: File
        inputBinding:
            prefix: "-V"
            position: 3
        secondaryFiles: [".tbi"]
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: "combined_somatic_plus_germline.vcf"
