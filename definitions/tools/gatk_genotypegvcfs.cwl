#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK HaplotypeCaller"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK-3.5.jar", "-T", "GenotypeGVCFs"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: mgibio/cle
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
