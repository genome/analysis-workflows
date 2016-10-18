#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add GT tags"
baseCommand: "/usr/bin/awk"
requirements:
    DockerRequirement:
        dockerPull: "ubuntu:xenial"
arguments:
    - '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}'
stdin: $(inputs.vcf.path)
stdout: 'output.gt.vcf'
inputs:
    vcf:
        type: File
outputs:
    processed_vcf:
        type: stdout

