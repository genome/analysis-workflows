#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "downsample unaligned BAM"

baseCommand: ["/gatk/gatk", "--java-options", "-Xmx16g", "DownsampleSam"]
arguments:
    ["--OUTPUT=", { valueFrom: $(runtime.outdir)/$(inputs.sam.nameroot).bam }, "--CREATE_INDEX", "--CREATE_MD5_FILE"]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.4.1"
inputs:
    sam:
        type: File
        inputBinding:
            prefix: "--INPUT="
            separate: false
    probability:
        type: float
        inputBinding:
            prefix: "--PROBABILITY="
            separate: false
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "--REFERENCE_SEQUENCE="
            separate: false
    random_seed:
        type: int?
        inputBinding:
            prefix: "--RANDOM_SEED="
            separate: false
    strategy:
        type: 
            - "null"
            - type: enum
              symbols: ["HighAccuracy", "ConstantMemory", "Chained"]
        inputBinding:
            prefix: "--STRATEGY="
            separate: false
outputs:
    downsampled_sam:
        type: File
        secondaryFiles: ['.md5', '^.bai']
        outputBinding:
            glob: $(inputs.sam.nameroot).bam
