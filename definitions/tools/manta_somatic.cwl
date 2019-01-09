#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Set up and execute manta using helper script"

requirements:
    - class: DockerRequirement
      dockerPull: mgibio/manta_somatic-cwl:1.5.0
    - class: InlineJavascriptRequirement
baseCommand: ["/usr/bin/python", "/usr/bin/manta_helper.py"]
arguments: 
    - position: 5
      valueFrom: $(runtime.outdir)
      prefix: "--runDir"
inputs:
    normal_bam:
        type: File?
        inputBinding:
            position: 2
            prefix: "--normalBam"
        secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
    tumor_bam:
        type: File
        inputBinding:
            position: 3
            prefix: "--tumorBam"
        secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
    reference:
        type: string
        inputBinding:
            position: 4
            prefix: "--referenceFasta"
    call_regions:
        type: File?
        inputBinding:
            position: 5
            prefix: "--callRegions"
    non_wgs:
        type: boolean?
        inputBinding:
            position: 6
            prefix: "--exome"
    output_contigs:
        type: boolean?
        inputBinding:
            position: 7
            prefix: "--outputContig"
outputs:
    diploid_variants:
        type: File
        outputBinding:
            glob: results/variants/diploidSV.vcf.gz
    somatic_variants:
        type: File
        outputBinding:
            glob: results/variants/somaticSV.vcf.gz
    all_candidates:
        type: File
        outputBinding:
            glob: results/variants/candidateSV.vcf.gz
    small_candidates:
        type: File
        outputBinding:
            glob: results/variants/candidateSmallIndels.vcf.gz
    tumor_only_variants:
        type: File?
        outputBinding:
            glob: results/variants/tumorSV.vcf.gz
