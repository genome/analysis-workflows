#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "som.py compares variants in vcf by location and alleles (using bcftools isec)."
baseCommand: ["/opt/hap.py/bin/som.py"]
requirements:
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "pkrusche/hap.py:v0.3.9"
arguments: [
    "-N",
    "-o", "sompy",
]
stdout: "sompy.out"
inputs:
    truth_vcf:
        type: File
        inputBinding:
            position: -2
    query_vcf:
        type: File
        inputBinding:
            position: -1
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-r"
            position: 1
    roi_bed:
        type: File
        doc: 'region of interest bed file to restrict'
        inputBinding:
            prefix: "-R"
            position: 2
outputs:
    sompy_out:
        type: stdout
