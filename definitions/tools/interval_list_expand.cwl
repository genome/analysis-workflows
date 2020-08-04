#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "expand interval list regions by a given number of basepairs"

baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "IntervalListTools"]

#/usr/bin/java -Xmx16g -jar /usr/picard/picard.jar IntervalListTools INPUT=94b23d66f667472690a7ee165e2037b6.interval_list OUTPUT=expanded.interval_list PADDING=100 UNIQUE=TRUE

arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/expanded.interval_list }, "UNIQUE=TRUE"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
inputs:
    interval_list:
        type: File
        inputBinding:
            prefix: "INPUT="
    roi_padding:
        type: int
        inputBinding:
            prefix: "PADDING="
outputs:
    expanded_interval_list:
        type: File
        outputBinding:
            glob: "expanded.interval_list"
