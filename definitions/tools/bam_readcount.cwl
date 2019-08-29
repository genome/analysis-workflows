#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bam-readcount"

baseCommand: ["/usr/bin/python", "/usr/bin/bam_readcount_helper.py"]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/bam_readcount_helper-cwl:1.1.1"
    - class: ResourceRequirement
      ramMin: 16000
    - class: InlineJavascriptRequirement
arguments: [
    { valueFrom: $(runtime.outdir), position: -3 }
]
stdout: $(inputs.sample)_bam_readcount.tsv
inputs:
    vcf:
        type: File
        inputBinding:
            position: -8
    sample:
        type: string
        inputBinding:
            position: -7
    reference_fasta:
        type: string
        inputBinding:
            position: -6
    bam:
        type: File
        inputBinding:
            position: -5
        secondaryFiles: [.bai]
    prefix:
        type: string?
        default: 'NOPREFIX'
        inputBinding:
            position: -4
    min_base_quality:
        type: int?
        default: 20
        inputBinding:
            position: -2
    min_mapping_quality:
        type: int?
        default: 0
        inputBinding:
            position: -1
outputs:
    snv_bam_readcount_tsv:
        type: File
        outputBinding:
            glob: |
                    ${
                        var name = "_bam_readcount_snv.tsv";
                        if (inputs.prefix.equals("NOPREFIX")) {
                            name = inputs.sample + name;
                        }
                        else {
                            name = inputs.prefix + "_" + inputs.sample + name;
                        }
                        return name;
                    }
    indel_bam_readcount_tsv:
        type: File
        outputBinding:
            glob: |
                    ${
                        var name = "_bam_readcount_indel.tsv";
                        if (inputs.prefix.equals("NOPREFIX")) {
                            name = inputs.sample + name;
                        }
                        else {
                            name = inputs.prefix + "_" + inputs.sample + name;
                        }
                        return name;
                    }
            
