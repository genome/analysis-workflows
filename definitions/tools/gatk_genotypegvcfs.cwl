#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK HaplotypeCaller"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "GenotypeGVCFs"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.5.0"
arguments:
    ["-o", 'genotype.vcf.gz']
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    gvcfs:
        type:
            type: array
            items: File
            inputBinding:
                prefix: "--variant"
        inputBinding:
            position: 2
outputs:
    genotype_vcf:
        type: File
        outputBinding:
            glob: "genotype.vcf.gz"
        secondaryFiles: [.tbi]
