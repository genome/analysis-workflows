#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants for Panel of Normals"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/cle
arguments:
    ["-minN", "2",
     "--setKey", "null",
     "--filteredAreUncalled",
     "--filteredrecordsmergetype", "KEEP_IF_ANY_UNFILTERED",
     "-o", { valueFrom: $(runtime.outdir)/panel_of_normals.vcf.gz }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    interval_list:
        type: File?
        inputBinding:
            prefix: "-L"
            position: 2
    normal_vcfs:
        type:
            type: array
            items: File
            inputBinding:
                prefix: "-V"
        secondaryFiles: [".tbi"]
        inputBinding:
            position: 3
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: "panel_of_normals.vcf.gz"

