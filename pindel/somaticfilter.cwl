#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel somatic filter v1"
baseCommand: ["/usr/bin/perl", "/usr/bin/somatic_indelfilter.pl", "filter.config"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: InitialWorkDirRequirement
      listing:
          - entryname: pindel_output_summary
            entry: $(inputs.pindel_output_summary)
          - entryname: 'filter.config'
            entry: |
                input=pindel_output_summary
                vaf=0.1
                cov=20
                hom=6
                pindel2vcf=/usr/bin/pindel2vcf
                reference=$(inputs.reference.path)
                referencename=GRCh38DH
                referencedate=20161216
                output=$(runtime.outdir)/pindel.out.vcf
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    pindel_output_summary:
        type: File
outputs:
    vcf:
        type: File
        outputBinding:
            glob: "pindel.out.vcf"
