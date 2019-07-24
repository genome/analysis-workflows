#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Sanitize a VCF"
baseCommand: ["/bin/bash","sanitize.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
      coresMin: 1
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'sanitize.sh'
        entry: |
            set -eou pipefail

            # remove lines containing non ACTGN bases, as they conflict with the VCF spec
            # and cause GATK to choke
            gunzip -c "$1" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' >sanitized.vcf
            /opt/htslib/bin/bgzip sanitized.vcf
            tabix sanitized.vcf.gz
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
outputs:
    sanitized_vcf:
        type: File
        outputBinding:
            glob: "sanitized.vcf.gz"
