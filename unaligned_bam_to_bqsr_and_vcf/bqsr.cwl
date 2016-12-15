#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'create BQSR table'
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator"]
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/bqsr.table },
    "--preserve_qscores_less_than", "6",
    "--disable_auto_index_creation_and_locking_when_reading_rods",
    "--disable_bam_indexing",
    "-dfrac", ".1",
    "-nct", "4"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/gatk-3.6:1"
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [.bai]
    known_sites:
        type:
            type: array
            items: File
            inputBinding:
                prefix: "-knownSites"
                position: 3
outputs:
    bqsr_table:
        type: File
        outputBinding:
            glob: "bqsr.table"
