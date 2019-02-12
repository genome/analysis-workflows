#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Ensembl Variant Effect Predictor"
baseCommand: ["/usr/bin/perl", "-I", "/opt/lib/perl/VEP/Plugins", "/usr/bin/variant_effect_predictor.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 64000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/vep_helper-cwl:1.0.0"
arguments:
    ["--format", "vcf",
    "--vcf",
    "--fork", "4",
    "--plugin", "Downstream",
    "--plugin", "Wildtype",
    "--term", "SO",
    "--transcript_version",
    "--offline",
    "--cache",
    "--everything",
    "-o", { valueFrom: $(runtime.outdir)/annotated.vcf }]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    cache_dir:
        type: string
        inputBinding:
            prefix: "--dir"
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
    pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
        default: "flag_pick"
        inputBinding:
            prefix: '--'
            separate: false
            position: 7
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
                        return []
                    }
                }
            position: 6
    custom_clinvar_vcf:
        type: File?
        secondaryFiles: [.tbi]
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.custom_clinvar_vcf) {
                        return ["--custom", inputs.custom_clinvar_vcf.path + ",clinvar,vcf,exact,0,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES"]

                    }
                    else {
                        return []
                    }
                }
            position: 7
    reference:
        type: string?
        inputBinding:
            prefix: "--fasta" 
            position: 8
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "annotated.vcf"
    vep_summary:
        type: File
        outputBinding:
            glob: "annotated.vcf_summary.html"
