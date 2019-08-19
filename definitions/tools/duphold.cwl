#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: "/usr/bin/duphold"
arguments: ["--threads", "$(runtime.cores)"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/duphold-cwl:0.1.4"

inputs:
    bam:
        type: File
        inputBinding:
            position: 1
            prefix: "--bam"
        doc: "aligned bam/cram file"
    output_vcf_name:
        type: string?
        default: "duphold_annotated.vcf"
        inputBinding:
            position: 4
            prefix: "--output"
        doc: "output vcf file name"
    reference:
        type: string
        inputBinding:
            position: 2
            prefix: "--fasta"
        doc: "reference used to align bam"
    snps_vcf:
        type: File?
        inputBinding:
            position: 3
            prefix: "--snp"
        doc: "snps vcf file to annotate hom/het variant counts within sv"
    sv_vcf:
        type: File
        inputBinding:
            position: 4
            prefix: "--vcf"
        doc: "sv vcf file to annotate"

outputs:
    annotated_sv_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

