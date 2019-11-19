#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Ensembl Variant Effect Predictor"
baseCommand: ["/usr/bin/perl", "-I", "/opt/lib/perl/VEP/Plugins", "/usr/bin/variant_effect_predictor.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/vep_custom_annotation.yml
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 64000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/vep_helper-cwl:1.1.0"
arguments:
    ["--format", "vcf",
    "--vcf",
    "--fork", "4",
    "--term", "SO",
    "--transcript_version",
    "--offline",
    "--cache",
    "--symbol",
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
    custom_annotations:
        type:
            - "null"
            - type: array
              items: ../types/vep_custom_annotation.yml#vep_custom_annotation
              label: "custom type, check types directory for input format"
              inputBinding:
                  valueFrom: |
                       ${
                           return [self.annotation.check_existing ? '--check_existing' : '',
                             '--custom',
                             [self.annotation.file.path,
                             self.annotation.name,
                             self.annotation.data_format,
                             self.method,
                             self.force_report_coordinates ? 1 : 0,
                             self.annotation.vcf_fields ? self.annotation.vcf_fields : ''
                             ].filter(String).join(',')
                           ].filter(String)
                       }
                  position: 6
    reference:
        type:
            - "null"
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .index]
        inputBinding:
            prefix: "--fasta" 
            position: 7
    plugins:
        type:
            type: array
            items: string
            inputBinding:
                prefix: "--plugin"
        inputBinding:
            position: 8
    everything:
        type: boolean?
        default: true
        inputBinding:
            prefix: "--everything"
            position: 9
    ensembl_assembly:
        type: string
        inputBinding:
            prefix: "--assembly"
            position: 10
        doc: "genome assembly to use in vep. Examples: 'GRCh38' or 'GRCm38'"
    ensembl_version:
        type: string
        inputBinding:
            prefix: "--cache_version"
            position: 11
        doc: "ensembl version - Must be present in the cache directory. Example: '95'"
    ensembl_species:
        type: string
        inputBinding:
            prefix: "--species"
            position: 12
        doc: "ensembl species - Must be present in the cache directory. Examples: 'homo_sapiens' or 'mus_musculus'"
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "annotated.vcf"
    vep_summary:
        type: File
        outputBinding:
            glob: "annotated.vcf_summary.html"
