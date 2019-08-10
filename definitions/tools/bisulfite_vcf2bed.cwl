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
            
            if [[ $3 == "no" ]]
            then
                #Creates a gzipped bed and a bedgraph that leaves out MT, random, GL contigs, etc
                /usr/bin/biscuit vcf2bed -t cg -k 1 -e $1 | /usr/bin/biscuit mergecg $2 /dev/stdin -k 2 |  tee >(/bin/gzip >cpgs.bed.gz) | cut -f 1-4 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpgs.bedgraph
            else
                /usr/bin/biscuit vcf2bed -t cg -k 1 -e $1 | /usr/bin/biscuit mergecg $2 /dev/stdin -k 2 |  tee >(/bin/gzip >cpgs.bed.gz) | cut -f 1-4 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpgs.bedgraph
                /usr/bin/biscuit vcf2bed -t ch -k 1 -e $1 | awk '$6 == "CA"' | tee >(/bin/gzip >cpas.bed.gz) | cut -f 1-3,8 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpas.bedgraph
                /usr/bin/biscuit vcf2bed -t ch -k 1 -e $1 | awk '$6 == "CT"' | tee >(/bin/gzip >cpts.bed.gz) | cut -f 1-3,8 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpts.bedgraph
                /usr/bin/biscuit vcf2bed -t ch -k 1 -e $1 | awk '$6 == "CC"' | tee >(/bin/gzip >cpcs.bed.gz) | cut -f 1-3,8 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpcs.bedgraph
            fi

inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    reference:
        type: string
        inputBinding:
            position: 2
    assay_non_cpg_sites:
        type: string
        inputBinding:
            position: 3
outputs:
    final_bed:
        type: File[]
        outputBinding:
            glob: "*.bed.gz"
    final_bedgraph:
        type: File[]
        outputBinding:
            glob: "*.bedgraph"
