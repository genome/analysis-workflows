#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "convert a biscuit pileup vcf to a bed and bedgraph file, via a wrapper script"
baseCommand: ["/bin/bash","bsvcf2bed.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 2
    - class: DockerRequirement
      dockerPull: "mgibio/biscuit:0.3.8"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'bsvcf2bed.sh'
        entry: |
            set -eou pipefail

            pileup_vcf="$1"
            reference="$2"
            cpgs_bed_out="$3"
            cpgs_bedgraph_out="$4"

            #Creates a gzipped bed and a bedgraph that leaves out MT, random, GL contigs, etc
            /usr/bin/biscuit vcf2bed -t cg -k 1 -e "$pileup_vcf" | /usr/bin/biscuit mergecg "$reference" /dev/stdin -k 2 |  tee >(/bin/gzip >"$cpgs_bed_out") | cut -f 1-4 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >"$cpgs_bedgraph_out"

arguments: [
    "$(runtime.outdir)/cpgs.bed.gz",
    "$(runtime.outdir)/cpgs.bedgraph"
]
inputs:
    vcf:
        type: File
        inputBinding:
            position: -2
    reference:
        type: string
        inputBinding:
            position: -1
outputs:
    cpgs:
        type: File
        outputBinding:
            glob: "cpgs.bed.gz"
    cpg_bedgraph:
        type: File
        outputBinding:
            glob: "cpgs.bedgraph"
