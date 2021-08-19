#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Sambamba: merge"
baseCommand: ["/usr/bin/perl", "merge.pl"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/bam-merge:0.1"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'merge.pl'
        entry: |
            #!/usr/bin/perl

            use strict;
            use warnings;

            use Getopt::Std;
            use File::Copy;

            my %opts;
            getopts('t:n:s', \%opts);
            my $nthreads = $opts{t} // die 'missing thread count';
            my $outfilename = $opts{n} // die 'missing output filename';
            my $sorted = $opts{s};

            my @bams = @ARGV;
            die 'missing input bams' unless scalar(@bams);

            #if there is only one bam, just copy it and index it
            if (scalar(@bams) == 1) {
                copy($bams[0], $outfilename) or die 'failed to copy file:' . $!;
            } else {
                if ($sorted) {
                    my $rv = system((qw(/usr/bin/sambamba merge -t)), $nthreads, $outfilename, @bams);
                    $rv == 0 or die 'failed to merge with sambamba';
                } else { #unsorted bams, use picard
                    my @args = (
                        'OUTPUT=' . $outfilename,
                        'ASSUME_SORTED=true',
                        'USE_THREADING=true',
                        'SORT_ORDER=unsorted',
                        'VALIDATION_STRINGENCY=LENIENT',
                        map { 'INPUT=' . $_ } @bams
                    );
                    my $rv = system((qw(java -jar -Xmx6g /opt/picard/picard.jar MergeSamFiles)), @args);
                    $rv == 0 or die 'failed to merge with picard';
                }
            }
            if ($sorted) {
                my $rv = system((qw(/usr/bin/sambamba index)), $outfilename);
                $rv == 0 or die 'failed to index';
            }

arguments: [
    "-t", "$(runtime.cores)",
]
inputs:
    bams:
        type: File[]
        inputBinding:
            position: 3
    sorted:
        type: boolean?
        default: false
        inputBinding:
            prefix: "-s"
            position: 2
    name:
        type: string?
        default: "merged.bam"
        inputBinding:
            prefix: "-n"
            position: 1
outputs:
    merged_bam:
        type: File
        outputBinding:
            glob: "$(inputs.name)"
