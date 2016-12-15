#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "mark duplicates"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "MarkDuplicates"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/DuplicatesMarked.bam },
    "ASSUME_SORT_ORDER=", "queryname",
    "METRICS_FILE=", "mark_dups_metrics.txt",
    "QUIET=", "true",
    "VALIDATION_STRINGENCY=", "LENIENT"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
            position: 1
outputs:
    duplicate_marked_bam:
        type: File
        outputBinding:
            glob: "DuplicatesMarked.bam"
