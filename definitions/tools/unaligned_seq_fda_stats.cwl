#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'unaligned_stats.pl']
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.16.1"
    - class: ResourceRequirement
      ramMin: 4000
    - class: StepInputExpressionRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'unaligned_stats.pl'
        entry: |
                #!/usr/bin/perl

                use strict;
                use warnings;
                use Carp;

                use FileHandle;
                use File::Basename;
                use File::Spec::Functions;
                use IO::Uncompress::AnyUncompress;


                # sets global variables with the default
                use vars qw/$filter $offset $samtools/;

                # defines the minimum base call quality score to filter
                # inclusive, for example, Bases with >= Q30
                $filter = 30;

                # defines the offset to calculate the base call quality score from the Phred +33 encoded by ASCII characters
                # For example, 33 = 63.06707094 - 30.06707094 (see below)
                $offset = 33;

                # specifies the program paths in docker
                $samtools = "/opt/samtools/bin/samtools";


                # main subroutine
                Main();

                # program exits here
                exit 0;


                sub Main {
                    # for the paths of input files
                    my @paths = @ARGV if @ARGV > 0;
                    croak "input file path required" unless @paths > 0;


                    # reads input files
                    my (%freq, %count, %base);
                    my $format = "fastq";
                    my $n = 0;
                    foreach my $path (@paths)
                    {
                        my $count_for_path = {}; # $count{$path};
                        my $freq_for_path = {};  # $freq{$path};

                        my $base_for_path = 0; # $base{$path};

                        my $fh;
                        if ($path =~ /\.gz$/ || $path =~ /\.bz2$/)
                        {
                            $fh = IO::Uncompress::AnyUncompress->new($path, MultiStream => 1);
                        }
                        elsif ($path =~ /\.bam$/)
                        {
                            # for MGI legacy instrument data
                            $format = "bam";

                            $fh = FileHandle->new("$samtools view $path |");
                        }
                        else
                        {
                            $fh = FileHandle->new($path, "r");
                        }


                        # reads an input file
                        croak "Cannot find a file: $path" unless -e $path;
                        croak "Cannot open a file: $path" unless defined $fh;
                        if ($format eq "fastq")
                        {
                            while (my $i = $fh->getline)
                            {
                                if ($n % 4 == 3)
                                {
                                    chomp $i;

                                    update_hash($i, $freq_for_path, $count_for_path);
                                }
                                elsif ($n % 4 == 1)
                                {
                                    chomp $i;

                                    my $nbase = count_N($i);

                                    $base_for_path += $nbase;
                                }


                                $n ++;
                            }
                        }
                        elsif ($format eq "bam")
                        {
                            while (my $i = $fh->getline)
                            {
                                next if index($i, '@') == 0;
                                chomp $i;

                                my @fields = split /\t/, $i;

                                update_hash($fields[10], $freq_for_path, $count_for_path);

                                my $nbase = count_N($fields[9]);

                                $base_for_path += $nbase;

                                $n ++;
                            }
                        }
                        else
                        {
                            croak "Invalid input file format: $format";
                        }

                        for my $freq (keys %$freq_for_path) {
                            $freq{sum}{$freq} += $freq_for_path->{$freq};
                        }
                        $freq{$path} = $freq_for_path;

                        for my $c (keys %$count_for_path) {
                            $count{sum}{$c} += $count_for_path->{$c};
                        }
                        $count{$path} = $count_for_path;

                        $base{sum} += $base_for_path;
                        $base{$path} = $base_for_path;

                        # closes the file handler
                        printf "\n%d lines (cumulative) read from %s", $n, $path;
                        $fh->close;
                    }


                    # prints out the source file information
                    printf "\n\n[Input file information: %d file(s)]", scalar(@paths);
                    print "\nfile\tdir";
                    foreach my $path (@paths)
                    {
                        my ($name, $dir) = fileparse($path);
                        printf "\n%s\t%s", $name, $dir;
                    }

                    foreach my $key ("sum", @paths)
                    {
                        # calculates the total number of base calls and sequences
                        my ($nseq, $ncall) = (0, 0);
                        foreach my $len (keys %{$count{$key}})
                        {
                            $ncall += $len * $count{$key}->{$len};
                            $nseq += $count{$key}->{$len};
                        }

                        croak "Exception: Invalid number of sequences, $nseq" unless $nseq == ($format eq "fastq") ? $n / 4 : $n;


                        # calculates the median and mean from a frequency hash
                        my ($median, $mean, $ncall2) = stat_freq($freq{$key});
                        croak "Exception: invalid count number of base calls" unless $ncall == $ncall2;

                        # counts by the minimum base call quality score
                        my $nfilter = 0;
                        foreach my $score (keys %{$freq{$key}})
                        {
                            if ($filter <= (ord($score) - $offset))
                            {
                                $nfilter += $freq{$key}->{$score};
                            }
                        }

                        # prints out the summary
                        printf "\n\n[Quality score summary from %s]", $key eq "sum" ? sprintf("%d file(s)", scalar(@paths)) : $key;
                        printf "\nTotal Number of Reads\t%s", $nseq;                                                                    # total sequence number
                        printf "\nTotal Number of basecalls\t%s", $ncall;
                        printf "\nBases with N\t%s\t%s (%%)", $base{$key}, $base{$key} / $ncall * 100;
                        printf "\nMedian Basecall Quality Score\t%s", $median;
                        printf "\nMean Basecall Quality Score\t%s", $mean;
                        printf "\nBases with >= Q%s\t%d\t%s (%%)", $filter, $nfilter, $nfilter / $ncall * 100;
                        print "\n";
                    }
                }


                sub update_hash {
                    my ($qual, $freq, $count) = @_;

                    # converts characters to ASCII numbers
                    my @scores = unpack('(a)*', $qual);

                    # updates the frequency of basecalling score
                    foreach (@scores)
                    {
                        $freq->{$_} ++;
                    }

                    # counts the sequence length
                    my $length = scalar(@scores);
                    $count->{$length} ++;
                }


                sub count_N {
                    my ($str) = @_;

                    # calculates the frequency of "N"
                    my $count = 0;
                    if (index($str, 'N') > -1)
                    {
                        $count = $str =~ tr/N//;
                    }

                    return $count;
                }

                sub stat_freq {
                    my ($hash) = @_;

                    # frequency hash: {key => data value, value => data count}
                    croak "empty frequency hash" unless keys(%$hash) > 0;

                    # sorts the data values
                    my @sort = sort {ord($a) <=> ord($b)} keys %$hash;

                    # counts the total number of data
                    my $n = 0;
                    map { $n += $hash->{$_} } @sort;

                    # calculates the median index (0-based index)
                    my @midx = ($n % 2 == 1) ? ((1 + $n) / 2) - 1 : ((1 + $n) / 2 - 1.5, (1 + $n) / 2 - 0.5);

                    # calculates the median of base call quality score from the frequency
                    my ($median, $mean);
                    my $idx = 0;
                    for (my $i=0; $i<@sort; $i ++)
                    {
                        # gets the last index of the given score in the array (1-based)
                        my $score_chr = $sort[$i];
                        my $score_num = ord($score_chr) - $offset;
                        $idx += $hash->{$score_chr};

                        # to calculate the mean
                        $mean += $score_num * $hash->{$score_chr};

                        # calculates the median
                        unless (defined $median)
                        {
                            if ($midx[0] <= $idx - 1)
                            {
                                if (@midx == 2)
                                {
                                    if ($midx[1] <= $idx - 1)
                                    {
                                        $median = $score_num;
                                    }
                                    else
                                    {
                                        $median = ($score_num + $sort[$i + 1]) / 2;
                                    }
                                }
                                else
                                {
                                    $median = $score_num;
                                }
                            }
                        }
                    }

                    return ($median, $mean / $n, $n);
                }
inputs:
    input_files:
        type: File[]
        inputBinding:
            position: 1
    output_name:
        type: string?
        default: $(inputs.input_files[0].basename)
stdout: "$(inputs.output_name)_unaligned_metrics.txt"
outputs:
    unaligned_stats:
        type: stdout
