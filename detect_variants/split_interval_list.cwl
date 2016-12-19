#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'helper.pl']
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
    - class: InitialWorkDirRequirement
      listing:
          - entryname: 'helper.pl'
            entry: |
                use File::Copy;

                system('/usr/bin/java', '-jar', '/usr/picard/picard.jar', 'IntervalListTools', @ARGV);

                my $i = 1;
                for(glob('*/scattered.interval_list')) {
                    #create unique names and relocate all the scattered intervals to a single directory
                    File::Copy::move($_, qq{$i.interval_list});
                    $i++
                }

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
