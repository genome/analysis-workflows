#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

arguments: ["/opt/AnnotSV_2.1/bin/AnnotSV", "-bedtools", "/usr/bin/bedtools", "-outputDir", "$(runtime.outdir)"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "mgibio/annotsv-cwl:2.1"

inputs:
    genome_build:
        type: string
        inputBinding:
            position: 2
            prefix: "-genomeBuild"
    input_vcf:
        type: File
        inputBinding:
            position: 3
            prefix: "-SVinputFile"
        doc: "vcf file to filter"
    output_tsv_name:
        type: string?
        default: "AnnotSV.tsv"
        inputBinding:
            position: 4
            prefix: "-outputFile"
        doc: "output file name"
    snps_vcf:
        type: File[]?
        inputBinding:
            position: 5
            prefix: "-vcfFiles"
            itemSeparator: ","
        doc: "snps vcf(s) for adding hom/het snp counts found within svs"

outputs:
    sv_variants_tsv:
        type: File
        outputBinding:
            glob: $(inputs.output_tsv_name)
