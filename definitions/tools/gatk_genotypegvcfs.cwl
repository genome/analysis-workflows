#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK GenotypeGVCFs"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx14g -Xms5g", "GenotypeGVCFs"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1" 
arguments:
    ["-G", "StandardAnnotation", "-O", 'genotype.vcf.gz']
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
    dbsnp_vcf:
        type: File?
        inputBinding:
            prefix: "--dbsnp"
            position: 3
        secondaryFiles: [.tbi]
outputs:
    genotype_vcf:
        type: File
        outputBinding:
            glob: "genotype.vcf.gz"
        secondaryFiles: [.tbi]
