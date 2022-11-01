#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'unaligned_stats.pl']
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.4.2"
    - class: ResourceRequirement
      ramMin: 16000
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


                # sets global variables with the default
                use vars qw/$filter $base $offset $bzcat $zcat $samtools/;

                # defines the minimum base call quality score to filter
                # inclusive, for example, Bases with >= Q30
                $filter = 30;

                # defines the basecalling letter to count
                # (default: "N")
                $base = 'N';

                # defines the offset to calculate the base call quality score from the Phred +33 encoded by ASCII characters
                # For example, 33 = 63.06707094 - 30.06707094 (see below)
                $offset = 33;

                # specifies the program paths in docker(mgibio/cle:v1.4.2)
                $bzcat = "/bin/bzcat";                                             # Version 1.0.6
                $zcat = "/bin/zcat";                                                # zcat (gzip) 1.6
                $samtools = "/opt/samtools/bin/samtools";              # samtools 1.3.1 using htslib 1.3.2


                # main subroutine
                Main();

                # program exits here
                exit 0;


                sub Main {
                    # gets the paths of input files
                    my @paths = @ARGV if @ARGV > 0;
                    croak "input file path required" unless @paths > 0;
                    
                    
                    # opens the input files
                    my (%freq, %count);
                    my $format = "fastq";
                    my $n = 0;                  # total number of lines
                    my $nbase = 0;
                    foreach my $path (@paths)
                    {
                        my $fh;
                        if ($path =~ /\.gz$/)
                        {
                            $fh = FileHandle->new("$zcat $path |");
                        }
                        elsif ($path =~ /\.bam$/)
                        {
                            $format = "bam";
                            
                            $fh = FileHandle->new("$samtools view $path |");
                        }
                        elsif ($path =~ /\.bz2$/)
                        {
                            $fh = FileHandle->new("$bzcat $path |");
                        }
                        else
                        {
                            # from a .sam or .fastq text file
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
                                    
                                    # updates the frequency/count hashes
                                    update_hash($i, \%freq, \%count, $offset);
                                }
                                elsif ($n % 4 == 1)
                                {
                                    chomp $i;
                                    
                                    # counts the frequency of a letter
                                    $nbase += count_base($i, $base);
                                }
                                
                                
                                $n ++;
                            }
                        }
                        elsif ($format eq "bam")
                        {
                            while (my $i = $fh->getline)
                            {
                                next if substr($i, 0, 1) eq '@';
                                chomp $i;
                                
                                my @fields = split /\t/, $i;
                                
                                # updates the frequency/count hashes
                                update_hash($fields[10], \%freq, \%count, $offset);
                                
                                # counts the frequency of a letter
                                $nbase += count_base($fields[9], $base);
                                
                                
                                $n ++;
                            }
                        }
                        else
                        {
                            croak "Invalid input file format: $format";
                        }
                        
                        # closes the file handler
                        printf "\n%d lines (cumulative) read from %s", $n, $path;
                        $fh->close;
                    }
                    
                    
                    # calculates the total number of base calls and sequences
                    my ($nseq, $ncall) = (0, 0);
                    foreach my $len (sort {$a <=> $b} keys %count)
                    {
                        $ncall += $len * $count{$len};
                        $nseq += $count{$len};
                    }
                    
                    croak "Exception: Invalid total number of sequences, $nseq" unless $nseq == ($format eq "fastq") ? $n / 4 : $n;
                    

                    # calculates the median and mean from a frequency hash
                    my ($median, $mean, $ncall2) = stat_freq(%freq);
                    croak "Exception: invalid count number of base calls" unless $ncall == $ncall2;
                    
                    # counts by the minimum base call quality score
                    my $nfilter = 0;
                    foreach my $score (keys %freq)
                    {
                        if ($filter <= $score)
                        {
                            $nfilter += $freq{$score};
                        }
                    }
                    
                    
                    # prints out the source file information
                    printf "\n\n[Input file information: %d file(s)]", scalar(@paths);
                    print "\nfile\tdir";
                    foreach my $path (@paths)
                    {
                        my ($name, $dir) = fileparse($path);
                        printf "\n%s\t%s", $name, $dir;
                    }
                    
                    # prints out the summary
                    printf "\n\n[Quality score summary]";
                    printf "\nTotal Number of Reads\t%s", $nseq;                                                                    # total sequence number
                    printf "\nTotal Number of basecalls\t%s", $ncall;
                    printf "\nBases with %s\t%s\t%s (%%)", $base, $nbase, $nbase / $ncall * 100;
                    printf "\nMedian Basecall Quality Score\t%s", $median;
                    printf "\nMean Basecall Quality Score\t%s", $mean;
                    printf "\nBases with >= Q%s\t%d\t%s (%%)", $filter, $nfilter, $nfilter / $ncall * 100;
                    print "\n\n";
                }


                sub update_hash {
                    my ($qual, $freq, $count, $offset) = @_;
                    
                    # converts characters to ASCII numbers
                    my @scores = map { unpack("C*", $_ ) - $offset } split("", $qual);
                    
                    # updates the hashes
                    map { $freq->{$_} ++ } @scores;
                    
                    # counts the sequence length
                    $count->{$#scores + 1} ++;
                }


                sub count_base {
                    my ($str, $base) = @_;
                    
                    # calculates the frequency of a letter
                    my $count = 0;
                    if (index($str, $base) > -1)
                    {
                        $count = $str =~ s/$base//g;
                    }
                    
                    return $count;
                }


                sub stat_freq {
                    my (%hash) = @_;
                    
                    # frequency hash: {key => data value, value => data count}
                    croak "empty frequency hash" unless keys(%hash) > 0;
                    
                    # sorts the data values
                    my @sort = sort {$a <=> $b} keys %hash;
                    
                    # counts the total number of data
                    my $n = 0;
                    map { $n += $hash{$_} } @sort;
                    
                    # calculates the median index (0-based index)
                    my @midx = ($n % 2 == 1) ? ((1 + $n) / 2) - 1 : ((1 + $n) / 2 - 1.5, (1 + $n) / 2 - 0.5);
                    
                    # calculates the median of base call quality score from the frequency
                    my ($median, $mean);
                    my $idx = 0;
                    for (my $i=0; $i<@sort; $i ++)
                    {
                        # gets the last index of the given score in the array (1-based)
                        my $score = $sort[$i];
                        $idx += $hash{$score};
                        
                        # to calculate the mean
                        $mean += $score * $hash{$score};
                        
                        # calculates the median
                        unless (defined $median)
                        {
                            if ($midx[0] <= $idx - 1)
                            {
                                if (@midx == 2)
                                {
                                    if ($midx[1] <= $idx - 1)
                                    {
                                        # same as ($score + $score) / 2
                                        $median = $score;
                                    }
                                    else
                                    {
                                        $median = ($score + $sort[$i + 1]) / 2;
                                    }
                                }
                                else
                                {
                                    $median = $score;
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
