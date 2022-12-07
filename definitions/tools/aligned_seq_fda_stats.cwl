#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'aligned_stats.pl']
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.4.2"
    - class: ResourceRequirement
      ramMin: 16000
    - class: StepInputExpressionRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'aligned_stats.pl'
        entry: |
                #!/usr/bin/perl
                
                use strict;
                use warnings;
                use Carp;
                
                use FileHandle;
                use File::Basename;
                use File::Spec::Functions; 
                
                
                # sets global variables with the default
                use vars qw/$samtools/;
                
                # specifies the program paths in docker(mgibio/cle:v1.4.2)
                $samtools = "/opt/samtools/bin/samtools";                # samtools 1.3.1 using htslib 1.3.2
                
                
                # main subroutine
                Main();
                
                # program exits here
                exit 0;
                
                
                sub Main {
                    # for the reference sequence FASTA file
                    my $pathfasta = shift @ARGV if @ARGV > 0;
                    croak "reference sequence FASTA required" unless defined $pathfasta && -e $pathfasta;
                    
                    # for the paths of input files
                    my @paths = @ARGV if @ARGV > 0;
                    croak "input file path required" unless @paths > 0;
                    
                    
                    # reads input files
                    my (%count);
                    my $n = 0;
                    foreach my $path (@paths)
                    {
                        croak "Invalid file path: $path" unless -e $path;
                        
                        my $fh;
                        if ($path =~ /\.bam$/)
                        {
                            $fh = FileHandle->new("$samtools view $path |");
                        }
                        elsif ($path =~ /\.cram$/)
                        {
                            $fh = FileHandle->new("$samtools view -T $pathfasta $path |");
                        }
                        else
                        {
                            $fh = FileHandle->new($path, "r");
                        }
                        
                        # reads a file
                        croak "Cannot open a file: $path" unless defined $fh;
                        while (my $i = $fh->getline)
                        {
                            next if index($i, '@') == 0;
                            chomp $i;
                            
                            my @fields = split /\t/, $i;
                            
                            # gets a FLAG value
                            my $flag = $fields[1];
                            
                            # parses the flag
                            my ($unmapped, $reverse, $first, $last, $secondary, $failed, $duplicate, $supplementary);
                            if ($flag & 0x4)
                            {
                               # "segment unmapped"
                               $unmapped = 1;
                            }
                            
                            if ($flag & 0x10)
                            {
                               # "being reverse complemented"
                               $reverse = 1;
                            }
                            
                            if ($flag & 0x40)
                            {
                                # "first in pair"
                                $first = 1;
                            }
                            
                            if ($flag & 0x80)
                            {
                                # "the last segment in the template"
                                $last = 1;
                            }
                            
                            if ($flag & 0x100)
                            {
                               # "secondary alignment"
                               $secondary = 1;
                            }
                            
                            if ($flag & 0x200)
                            {
                               # "not passing quality controls"
                               $failed = 1;
                            }
                            
                            if ($flag & 0x400)
                            {
                               # "PCR or optical duplicate"
                               $duplicate = 1;
                            }
                            
                            if ($flag & 0x800)
                            {
                               # "supplementary alignment"
                               $supplementary = 1;
                            }
                            
                            
                            # sorts reads by the flag information
                            if ($failed)
                            {
                                $count{failed} ++;
                            }
                            elsif ($duplicate)
                            {
                                $count{duplicate} ++;
                                
                                if ($unmapped)
                                {
                                    $count{"duplicate\tunmapped"} ++;
                                }
                                else
                                {
                                    if ($secondary || $supplementary)
                                    {
                                        $count{"duplicate\tfiltered"} ++;
                                    }
                                    else
                                    {
                                        $count{"duplicate\tprimary"} ++;
                                        
                                        if ($reverse)
                                        {
                                            $count{"duplicate\tprimary\treverse"} ++;
                                        }
                                        else
                                        {
                                            $count{"duplicate\tprimary\tforward"} ++;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                $count{unique} ++;
                                
                                if ($unmapped)
                                {
                                    $count{"unique\tunmapped"} ++;
                                    
                                    if ($first)
                                    {
                                        $count{"unique\tunmapped\tfirst"} ++;
                                    }
                                    elsif ($last)
                                    {
                                        $count{"unique\tunmapped\tlast"} ++;
                                    }
                                    else
                                    {
                                        croak "Exception in first and last: $i"
                                    }
                                }
                                else
                                {
                                    if ($secondary || $supplementary)
                                    {
                                        $count{"unique\tfiltered"} ++;
                                    }
                                    else
                                    {
                                        $count{"unique\tprimary"} ++;
                                        
                                        if ($reverse)
                                        {
                                            $count{"unique\tprimary\treverse"} ++;
                                        }
                                        else
                                        {
                                            $count{"unique\tprimary\tforward"} ++;
                                        }
                                        
                                        if ($first)
                                        {
                                            $count{"unique\tprimary\tfirst"} ++;
                                        }
                                        elsif ($last)
                                        {
                                            $count{"unique\tprimary\tlast"} ++;
                                        }
                                        else
                                        {
                                            croak "Exception in first and last: $i"
                                        }
                                    }
                                }
                            }
                            
                            unless ($secondary || $supplementary)
                            {
                                # total sequencing read number matching with that from raw fastq files if there are no reads filtered out by the aligner
                                $n ++;
                            }
                        }
                        
                        # closes the file handler
                        printf "\n%d reads (cumulative) from %s", $n, $path;
                        $fh->close;
                    }
                    
                    
                    # for missing values
                    my $default_zero = 0;
                    $count{"duplicate\tunmapped"} = $default_zero unless exists $count{"duplicate\tunmapped"};
                    $count{failed} = $default_zero unless exists $count{failed};
                    
                    
                    # prints out the source file information
                    printf "\n\n[Input file information: %d file(s)]", scalar(@paths);
                    print "\nfile\tdir";
                    foreach my $path (@paths)
                    {
                        my ($name, $dir, $ext) = fileparse($path, qr/\.[^.]*/);
                        printf "\n%s\t%s", $name . $ext, $dir;
                    }
                    
                    # prints out the summary
                    printf "\n\n[Flag summary from %d file(s)]", scalar(@paths);
                    printf "\nTotal Read Count (R1 + R2)\t%s", $n;
                    printf "\nQC-failed Read Count\t%s", $count{failed};
                    printf "\nUnique Read Pairs\t%s\t%s (%%)", $count{"unique\tprimary\tfirst"} + $count{"unique\tunmapped\tfirst"}, ($count{"unique\tprimary\tfirst"} + $count{"unique\tunmapped\tfirst"}) / $n * 100 * 2;
                    printf "\nTotal Mapped Reads\t%s\t%s (%%)", $count{"duplicate\tprimary"} + $count{"unique\tprimary"}, ($count{"duplicate\tprimary"} + $count{"unique\tprimary"}) / $n * 100;
                    printf "\nNon-Mapped Reads\t%s\t%s (%%)", $count{"duplicate\tunmapped"} + $count{"unique\tunmapped"}, ($count{"duplicate\tunmapped"} + $count{"unique\tunmapped"}) / $n * 100;
                    printf "\nUnique Mapped Reads\t%s\t%s (%%)", $count{"unique\tprimary"}, $count{"unique\tprimary"} / $n * 100;
                    printf "\nMapped Read Duplication\t%s\t%s (%%)", $count{"duplicate\tprimary"}, $count{"duplicate\tprimary"} / $n * 100;
                    printf "\nStrand ratio (forward, reverse, reverse/forward of unique mapped)\t%s\t%s\t%s", $count{"unique\tprimary\tforward"}, $count{"unique\tprimary\treverse"}, $count{"unique\tprimary\treverse"} / $count{"unique\tprimary\tforward"};
                    print "\n\n";
                }

inputs:
    reference:
        type:
            - string
            - File
        inputBinding:
            position: 1
    input_files:
        type: File[]
        inputBinding:
            position: 2
    output_name:
        type: string?
        default: $(inputs.input_files[0].basename)
stdout: "$(inputs.output_name)_aligned_metrics.txt"
outputs:
    aligned_stats:
        type: stdout

