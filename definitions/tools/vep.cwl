#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Ensembl Variant Effect Predictor"
baseCommand: ["/usr/bin/perl", "-I", "/opt/lib/perl/VEP/Plugins", "/usr/bin/variant_effect_predictor.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments:
    ["--format", "vcf",
    "--vcf",
    "--plugin", "Downstream",
    "--plugin", "Wildtype",
    "--symbol",
    "--term", "SO",
    "--flag_pick",
    "--transcript_version",
    "-o", { valueFrom: $(runtime.outdir)/annotated.vcf }]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    cache_dir:
        type: string?
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.cache_dir) {
                        return ["--offline", "--cache", "--dir", inputs.cache_dir ]
                    }
                    else {
                        return "--database"
                    }
                }
            position: 4
    synonyms_file:
        type: File?
        inputBinding:
            prefix: "--synonyms"
            position: 2
    coding_only:
        type: boolean
        inputBinding:
            prefix: "--coding_only"
            position: 3
        default: false
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.custom_gnomad_vcf) {
                        return ['--check_existing', '--custom', inputs.custom_gnomad_vcf.path + ',gnomADe,vcf,exact,0,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS']
                    }
                    else {
                        if (inputs.cache_dir) {
                            return ['--max_af', '--af_gnomad', '--af_1kg']
                        }
                        else {
                            return []
                        }
                    }
                }
            position: 6
    hgvs:
        type: boolean?
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.hgvs) {
                        if (inputs.cache_dir) {
                            return ["--hgvs", "--fasta", inputs.reference]
                        }
                        else {
                            return ["--hgvs"]
                        }
                    }
                    else {
                        return []
                    }
                }
            position: 5
    reference:
        type: string?
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "annotated.vcf"
    vep_summary:
        type: File
        outputBinding:
            glob: "annotated.vcf_summary.html"
