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

            die "wrong number of inputs" unless scalar(@ARGV) == 3;
            my ($output_dir, $interval_list, $scatter_count) = @ARGV;
            
            my $i = 1;

            if ($scatter_count == 1) {
                File::Copy::copy($interval_list,qq{$i.interval_list});
            } else {

                my $retval = system('/usr/bin/java', '-jar', '/usr/picard/picard.jar', 'IntervalListTools', 'OUTPUT='.$output_dir, 'INPUT='.$interval_list, 'SCATTER_COUNT='. $scatter_count);
                exit $retval if $retval != 0;

                for (glob('*/scattered.interval_list')) {
                    #create unique names and relocate all the scattered intervals to a single directory
                    File::Copy::move($_, qq{$i.interval_list});
                    $i++
                }
            }

arguments:
    [{ valueFrom: $(runtime.outdir) }]
inputs:
    interval_list:
        type: File
        inputBinding:
            position: 1
    scatter_count:
        type: int
        inputBinding:
            position: 2
outputs:
    split_interval_lists:
        type: File[]
        outputBinding:
            glob: "*.interval_list"
