#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "GATK HaplotypeCaller"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK-3.5.jar", "-T", "HaplotypeCaller"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: mgibio/cle
arguments:
    ["-o", { valueFrom: '${
            if (inputs.intervals.length == 1 && inputs.intervals[0].match(/^[0-9A-Za-z]+$/)) {
                return inputs.intervals[0] + ".g.vcf.gz";
            } else {
                return "output.g.vcf.gz";
            }
        }'
    }]
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
outputs:
    gvcf:
        type: File
        outputBinding:
            glob: "*.g.vcf.gz"
        secondaryFiles: [.tbi]
