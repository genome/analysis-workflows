#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'split_interval_list_helper.pl']
requirements:
    - class: ResourceRequirement
      ramMin: 6000
    - class: DockerRequirement
      dockerPull: broadinstitute/picard:2.24.2
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'split_interval_list_helper.pl'
        entry: |
            use File::Copy;

            my $i = 1;

            my ($scatter_count_arg) = grep{ /^SCATTER_COUNT=\d+$/ } @ARGV;

            # If only 1 scatter count is defined, only copy the file
            if ($scatter_count_arg && $scatter_count_arg =~ /^SCATTER_COUNT=(\d+)$/) {
                my $scatter_count = $1;
                if ($scatter_count == 1) {
                    my ($input_arg) = grep{ /^INPUT=.*$/ } @ARGV;
                    if ($input_arg && $input_arg =~ /^INPUT=(.*)$/) {
                        my $input = $1;
                        File::Copy::move($input,qq{$i.interval_list});
                    }
                }
            } else {

                my $retval = system('/usr/bin/java', '-jar', '/usr/picard/picard.jar', 'IntervalListTools', @ARGV);
                exit $retval if $retval != 0;

                for(glob('*/scattered.interval_list')) {
                    #create unique names and relocate all the scattered intervals to a single directory
                    File::Copy::move($_, qq{$i.interval_list});
                    $i++
                }
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
