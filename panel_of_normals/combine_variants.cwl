#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
baseCommand: ["/usr/local/bin/jdk1.8.0_45/bin/java", "-jar", "/usr/local/bin/GATK3.6/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: DockerRequirement
      dockerPull: "dbmi/gatk-docker:v1" #GATK 3.6 at a specific container revision
arguments:
    ["-setKey", "null",
     "--filteredAreUncalled", 
     "--filteredrecordsmergetype", "KEEP_IF_ANY_UNFILTERED",
     "-o", { valueFrom: $(runtime.outdir)/panel_of_normals.vcf.gz }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
        secondaryFiles: [".fai", "^.dict"]
    vcfs:
        type: 
            type: array
            items: File
        inputBinding:
            prefix: "--variant"
            position: 2
        secondaryFiles: [.tbi]
    minimumN:
        type: int
        inputBinding:
            prefix: "--minimumN"
            position: 3
    interval_list:
        type: File?
        inputBinding:
            prefix: "-L"
            position: 4
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: "panel_of_normals.vcf.gz"
        secondaryFiles: [.tbi]

