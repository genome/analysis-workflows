#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK HaplotypeCaller"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
requirements:
    - class: ResourceRequirement
      ramMin: 10000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.5.0"
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    cram:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [^.crai]
    emit_reference_confidence:
        type: string
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
    output_file_name:
        type: string
        default: "output.g.vcf.gz"
        inputBinding:
            prefix: "-o"
            position: 8
outputs:
    gvcf:
        type: File
        outputBinding:
            glob: $(inputs.output_file_name)
        secondaryFiles: [.tbi]
