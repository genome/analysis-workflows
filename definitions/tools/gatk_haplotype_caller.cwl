#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK HaplotypeCaller"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx16g", "HaplotypeCaller"]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [^.bai]
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
        inputBinding:
            prefix: "-ERC"
            position: 3
    gvcf_gq_bands:
        type:
            type: array
            items: string
            inputBinding:
                prefix: "-GQB"
        inputBinding:
            position: 4
    intervals:
        type:
            type: array
            items: string
            inputBinding:
                prefix: "-L"
        inputBinding:
            position: 5
    dbsnp_vcf:
        type: File?
        inputBinding:
            prefix: "--dbsnp"
            position: 6
        secondaryFiles: [.tbi]
    contamination_fraction:
        type: string?
        inputBinding:
            prefix: "-contamination"
            position: 7
    max_alternate_alleles:
        type: int?
        doc: 'maximum number of alternate alleles to genotype'
        inputBinding:
            prefix: '--max_alternate_alleles'
            position: 8
    ploidy:
        type: int?
        doc: 'number of chromosomes per sample'
        inputBinding:
            prefix: '-ploidy'
            position: 9
    read_filter:
        type: string?
        doc: 'filters to apply to reads before analysis'
        inputBinding:
            prefix: '--read_filter'
            position: 10
    output_file_name:
        type: string
        default: "output.g.vcf.gz"
        inputBinding:
            prefix: "-O"
            position: 11
outputs:
    gvcf:
        type: File
        outputBinding:
            glob: $(inputs.output_file_name)
        secondaryFiles: [.tbi]
