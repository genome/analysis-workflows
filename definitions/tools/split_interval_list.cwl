#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'split_interval_list_helper.pl']
requirements:
    - class: ResourceRequirement
      ramMin: 6000
    - class: DockerRequirement
      dockerPull: mgibio/cle:v1.3.1
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'split_interval_list_helper.pl'
        entry: |
            use File::Copy;

            my $retval = system('/usr/bin/java', '-jar', '/usr/picard/picard.jar', 'IntervalListTools', @ARGV);
            exit $retval if $retval != 0;

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
