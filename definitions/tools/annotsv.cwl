#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

arguments: ["/opt/AnnotSV_2.3/bin/AnnotSV", "-bedtools", "/usr/bin/bedtools", "-outputDir", "$(runtime.outdir)",  "-outputFile", "$(inputs.output_base).tsv"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "mgibio/annotsv-cwl:2.3"

inputs:
    genome_build:
        type: string
        inputBinding:
            position: 2
            prefix: "-genomeBuild"
        doc: "genome build used, GRCh37(tool default), GRCh38, mm9, or mm10"
    input_vcf:
        type: File
        inputBinding:
            position: 3
            prefix: "-SVinputFile"
        doc: "vcf file to filter"
    output_base:
        type: string?
        default: "AnnotSV"
        inputBinding:
        doc: "base for output file name"
    snps_vcf:
        type: File[]?
        inputBinding:
            position: 5
            prefix: "-snvIndelFiles"
            itemSeparator: ","
        doc: "snps vcf(s) for adding hom/het snp counts found within svs"
    annotations:
        type:
            - string
            - Directory
        inputBinding:
            position: 6
            prefix: "-annotationsDir"
        doc: "directory/path of the annotsv annotations directory"
outputs:
    annotated_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.output_base).tsv"
    unannotated_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.output_base).unannotated.tsv"
