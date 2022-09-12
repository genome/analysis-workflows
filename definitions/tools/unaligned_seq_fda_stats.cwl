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

                # It is free software; you can redistribute it and/or modify it under the terms of
                # the GNU General Public License (GPLv3) as published by the Free Software Foundation.
                #
                # McDonnell Genome Institute
                # Washington University School of Medicine
                # 4444 Forest Park Ave
                # St. Louis, Missouri 63108
                # 
                # http://genome.wustl.edu


                # includes perl modules
                use strict;
                use warnings;
                use Carp;

                use FileHandle;
                use File::Basename;
                use File::Spec::Functions; 
                use Cwd;


                # predeclares global variable names
                use vars qw/$version $about $is_windows $scriptname $scriptdir/;

                # finds the perl script name and path
                $scriptname = basename($0, '.pl') . '.pl';
                $scriptdir = getcwd;       # dirname($0);



                ## ------------ About the program ------------------------------------- #

                # the version
                $version = '2.0';

                $about = qq!
                FASTQ statistics $version
                Usage: $scriptname [FORMAT] [FILE1] [FILE2] [FILE3] ...
                Print out statistics calculated from FILE(s).
                FORMAT is required, either "fastq" or "bam" for FILE(s).

                Example:
                  \$ $scriptname fastq read1.fastq.gz read2.fastq.gz
                !;



                ## ------------ Setting global variables ------------------------------------- #

                # predeclares global variable names
                use vars qw/$format @paths $filter $bzcat $zcat $samtools $base $offset/;

                # defines the format of input files
                # (required, "fastq" or "bam")
                $format = "fastq";
                #$format = "bam";

                # lists the full paths of fastq or BAM files
                @paths = (
                    # H_NJ-HCC1395-HCC1395_BL
                    #"/gscmnt/gc2560/core/instrument_data/2895499331/csf_150397221/gerald_H7HY2CCXX_3_TGACCACG.bam",
                    #"/gscmnt/gc2560/core/instrument_data/2895499399/csf_150397349/gerald_H7HY2CCXX_4_TGACCACG.bam",
                );

                # defines the minimum base call quality score to filter
                # inclusive, for example, Bases with >= Q30
                $filter = 30;


                # specifies the program paths
                # docker(mgibio/cle:v1.4.2)
                $bzcat = "/bin/bzcat";                                             # Version 1.0.6
                $zcat = "/bin/zcat";                                                # zcat (gzip) 1.6
                $samtools = "/opt/samtools/bin/samtools";              # samtools 1.3.1 using htslib 1.3.2

                # specifies the N letter to count
                # (default: "N")
                $base = 'N';

                # defines the offset to calculate the base call quality score from the Phred +33 encoded by ASCII characters
                # For example, 33 = 63.06707094 - 30.06707094 (see below)
                $offset = 33;



                ## ------------ Main subroutine ---------------------------------------------- #
                # prints the program information
                print($about);

                # main subroutine
                Main();

                # this program is terminated here
                exit 0;



                ## ------------ Library of the main subroutine ----------------------------- #

                sub Main {
                    # local variables
                    
                    
                    # gets the input file format
                    unless (defined $format)
                    {
                        if (@ARGV > 0)
                        {
                            $format = shift @ARGV;
                        }
                        else
                        {
                            croak "FORMAT required: fastq or bam"
                        }
                        
                        croak "Invalid input file format: $format" unless $format eq "fastq" || $format eq "bam";
                    }
                    
                    # gets the full input paths
                    unless (@paths > 0)
                    {
                        @paths = @ARGV if (@ARGV > 0);
                        
                        croak "paths required" unless @paths > 0;
                    }
                    
                    
                    # opens input files
                    my $n = 0;                  # total number of lines
                    my $nbase = 0;
                    my (%freq, %count);
                    foreach my $path (@paths)
                    {
                        # creates a file handler
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
                        
                        
                        # reads a file
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
                                
                                # gets a QUAL value
                                # my $qual = $fields[10];
                                
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
                        
                        
                        #
                        printf "\n  %d lines (cumulative) read from %s", $n, $path;
                        
                        
                        # closes the file handler
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
                    
                    # filters with the minimum base call quality score
                    my $nfilter = 0;
                    foreach my $score (keys %freq)
                    {
                        if ($filter <= $score)
                        {
                            $nfilter += $freq{$score};
                        }
                    }
                    
                    
                    # prints out the source file information
                    printf "\n\n[Input file information]";
                    foreach my $path (@paths)
                    {
                        my ($name, $dir, $ext) = fileparse($path, qr/\.[^.]*/);
                        printf "\n%s\t%s", $name . $ext, $dir;
                    }
                    
                    
                    # prints out the summary
                    printf "\n\n[Quality score summary]";
                    printf "\nTotal Number of Reads\t%s", $nseq;                                                                    # total sequence number
                    printf "\nTotal Number of basecalls\t%s", $ncall;
                    printf "\nBases with %s\t%s\t%s (%%)", $base, $nbase, $nbase / $ncall * 100;
                    printf "\nMedian Basecall Quality Score\t%s", $median;
                    printf "\nMean Basecall Quality Score\t%s", $mean;
                    printf "\nBases with >= Q%s\t%d\t%s (%%)", $filter, $nfilter, $nfilter / $ncall * 100;
                    
                    
                    # Base call score frequency
                    printf "\n\n[Base call score frequency: %d sequences from %d file(s)]\n", $nseq, scalar(@paths);
                    printhash(%freq);
                    
                    # Sequence length statistics
                    printf "\n\n[Sequence length statistics]\n";
                    printhash(%count);
                    
                    
                    # the end of the main subroutine
                }



                ## ------------ Library of subroutines --------------------------------------- #


                sub update_hash {
                    # updates the frequency/count hashes
                    my ($qual, $freq, $count, $offset) = @_;

                    # local variables
                    
                    
                    # converts characters to ASCII numbers
                    my @scores = map { unpack("C*", $_ ) - $offset } split("", $qual);
                    
                    
                    # updates the hashes
                    map { $freq->{$_} ++ } @scores;
                    
                    #map { $count{$_} ++ }  split("", $i);
                    $count->{$#scores + 1} ++;                 # for the sequence length
                }


                sub count_base {
                    # counts the frequency of a letter
                    my ($str, $base) = @_;
                    
                    # local variables
                    
                    
                    # calculates the frequency of a letter
                    my $count = 0;
                    if (index($str, $base) > -1)
                    {
                        $count = $str =~ s/$base//g;
                    }
                    
                    
                    return $count;
                }


                sub stat_freq {
                    # calculates the median and mean from a frequency hash
                    my (%hash) = @_;
                    
                    # local variables
                    
                    # as default
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
                        
                        # calculates the mean and mode
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


                sub printhash {
                    # prints a hash
                    my (%hash) = @_;

                    # local variables
                    
                    foreach my $key (sort {$a <=> $b} keys %hash)
                    {
                        printf "%s\t%s\n", $key, $hash{$key};
                    }
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
