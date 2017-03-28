#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "mutect2 (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "MuTect2"]
requirements:
    - class: ResourceRequirement
      ramMin: 20000
      tmpdirMin: 100000
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/output.vcf.gz }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
        secondaryFiles: [".fai", "^.dict"]
    tumor_cram:
        type: File
        inputBinding:
            prefix: "-I:tumor"
            position: 2
        secondaryFiles: [^.crai]
    normal_cram:
        type: File?
        inputBinding:
            prefix: "-I:normal"
            position: 3
        secondaryFiles: [^.crai]
    interval_list:
        type: File
        inputBinding:
            prefix: "-L"
            position: 4
    dbsnp_vcf:
        type: File?
        inputBinding:
            prefix: "--dbsnp"
            position: 5
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        inputBinding:
            prefix: "--cosmic"
            position: 6
        secondaryFiles: [.tbi]
    artifact_detection_mode:
        type: boolean?
        inputBinding:
            prefix: "--artifact_detection_mode"
            position: 7 
outputs:
    vcf:
        type: File
        outputBinding:
            glob: "output.vcf.gz"
        secondaryFiles: [.tbi]
