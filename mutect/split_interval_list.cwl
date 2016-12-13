#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', '-MFile::Copy', '-e',
"system('/usr/bin/java', '-jar', '/usr/picard/picard.jar', 'IntervalListTools', @ARGV); my $i = 1; for(glob('*/scattered.interval_list')) { File::Copy::move($_, qq{$i.interval_list}); $i++ }"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
arguments:
    [{ valueFrom: OUTPUT=$(runtime.outdir) }]
inputs:
    interval_list:
        type: File
        inputBinding:
            prefix: "INPUT="
            separate: false
            position: 1
    scatter_count:
        type: int
        inputBinding:
            prefix: "SCATTER_COUNT="
            separate: false
            position: 2
outputs:
    split_interval_lists:
        type: File[]
        outputBinding:
            glob: "*.interval_list"
