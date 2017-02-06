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
    "-L", "chr1",
    "-L", "chr2",
    "-L", "chr3",
    "-L", "chr4",
    "-L", "chr5",
    "-L", "chr6",
    "-L", "chr7",
    "-L", "chr8",
    "-L", "chr9",
    "-L", "chr10",
    "-L", "chr11",
    "-L", "chr12",
    "-L", "chr13",
    "-L", "chr14",
    "-L", "chr15",
    "-L", "chr16",
    "-L", "chr17",
    "-L", "chr18",
    "-L", "chr19",
    "-L", "chr20",
    "-L", "chr21",
    "-L", "chr22",
    "-nct", "4"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/gatk-3.6:1"
    - class: ResourceRequirement
      ramMin: 16000
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
