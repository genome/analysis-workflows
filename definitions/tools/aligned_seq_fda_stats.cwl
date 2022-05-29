#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', 'aligned_stats.pl']
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-52"
    - class: ResourceRequirement
      ramMin: 16000
    - class: StepInputExpressionRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'aligned_stats.pl'
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
		$version = '1.6';

		$about = qq!
		BAM flag statistics $version
		Usage: $scriptname [FASTA] [FILE1] [FILE2] [FILE3] ...
		Print out flag statistics calculated from FILE(s).
		FASTA is the reference sequence fastq file.
		FILE(s) format can be SAM, BAM, and CRAM.

		Example:
		  \$ $scriptname all_sequences.fa normal.bam tumor.cram
		!;



		## ------------ Setting global variables ------------------------------------- #

		# predeclares global variable names
		use vars qw/@paths $samtools $pathfasta/;

		# lists the full paths of fastq or BAM files
		# (default: empty from @ARGV)
		@paths = (
		    # H_NJ-HCC1395-HCC1395_BL
		    #'/storage1/fs1/bga/Active/johnegarza/fastqc_hcc1395/DNA_normal/aligned/normal.bam',
		    #'/storage1/fs1/bga/Active/shared/gmsroot/model_data/80beed84e3104595862e7a2c7b7f896e/build8b9e77b87b0144b8a025288e92434686/tmp/cromwell-executions/immuno.cwl/58c62cbc-026b-4c9b-8d21-04b0aa4a20e8/call-somatic/somatic_exome.cwl/b24174a4-a185-4fe4-a222-50cfb2ffb305/call-normal_index_cram/inputs/1389193938/normal.bam.cram',
		    
		    #'/storage1/fs1/bga/Active/shared/gmsroot/model_data/80beed84e3104595862e7a2c7b7f896e/build191d8f37d22a4d79b2591c01cb80948e/results/normal.bam.cram'
		);


		# specifies the reference sequence FASTA file
		# (default: undef from @ARGV)
		$pathfasta = undef;
		#$pathfasta = '/gscmnt/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa';

		# specifies the program paths
		# docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-52)
		$samtools = "/usr/bin/samtools";                                  # Version: 1.10 (using htslib 1.10.2-3)



		## ------------ Main subroutine ---------------------------------------------- #
		# prints the program information
		print($about);

		# main subroutine
		Main();

		# this program is terminated here
		exit 1;



		## ------------ Library of the main subroutine ----------------------------- #

		sub Main {
		    # local variables
		    my $default_zero = 0;
		    
		    
		    # gets the input file format
		    if (defined $pathfasta)
		    {
			croak "Invalid reference sequence path: $pathfasta" unless -e $pathfasta;
		    }
		    else
		    {
			$pathfasta = shift @ARGV if (@ARGV > 0);
			
			croak "Invalid reference sequence path: $pathfasta" unless -e $pathfasta;
		    }
		    
		    # gets the full input paths
		    unless (@paths > 0)
		    {
			@paths = @ARGV if (@ARGV > 0);
			
			croak "paths required" unless @paths > 0;
		    }
		    
		    
		    #
		    my $n = 0;                  # total number of lines
		    my (%count);
		    foreach my $path (@paths)
		    {
			# creates a file handler
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
			    # from a .sam or .fastq text file
			    $fh = FileHandle->new($path, "r");
			}
			
			
			# reads a file
			croak "Cannot find a file: $path" unless -e $path;
			croak "Cannot open a file: $path" unless defined $fh;
			while (my $i = $fh->getline)
			{
			    next if substr($i, 0, 1) eq '@';
			    
			    chomp $i;
			    
			    my @fields = split /\t/, $i;
			    
			    # gets a FLAG value
			    #my $qname = $fields[0];
			    my $flag = $fields[1];
			    
			    # parses the flag
			    my ($unmapped, $reverse, $first, $last, $secondary, $failed, $duplicate, $supplementary);
			    if ($flag & 0x4)
			    {
			       # print "segment unmapped\n";
			       $unmapped = 1;
			    }
			    
			    if ($flag & 0x10)
			    {
			       #print "SEQ being reverse complemented\n";
			       $reverse = 1;
			    }
			    
			    if ($flag & 0x40)
			    {
				#print "first in pair\n";
				$first = 1;
			    }
			    
			    if ($flag & 0x80)
			    {
				#print "the last segment in the template\n";
				$last = 1;
			    }
			    
			    if ($flag & 0x100)
			    {
			       #print "secondary alignment\n";
			       $secondary = 1;
			    }
			    
			    if ($flag & 0x200)
			    {
			       #print "not passing quality controls\n";
			       $failed = 1;
			    }
			    
			    if ($flag & 0x400)
			    {
			       #print "PCR or optical duplicate\n";
			       $duplicate = 1;
			    }
			    
			    if ($flag & 0x800)
			    {
			       #print "supplementary alignment\n";
			       $supplementary = 1;
			    }
			    
			    
			    # sorts reads
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
				    # counts only the primary non-duplicate alignment 
				    if ($secondary || $supplementary)
				    {
					# skipped
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
				    # counts only the primary non-duplicate alignment 
				    if ($secondary || $supplementary)
				    {
					# skipped
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
				$n ++;        # total sequencing read number matching that from fastq files
			    }
			}
			
			
			# prints out a progress message
			printf "\n  %d reads (cumulative) from %s", $n, $path;
			
			
			# closes the file handler
			$fh->close;
		    }
		    
		    
		    # for missing values
		    $count{"duplicate\tunmapped"} = $default_zero unless exists $count{"duplicate\tunmapped"};
		    $count{failed} = $default_zero unless exists $count{failed};
		    
		    
		    # prints out the summary
		    printf "\n\n[Flag summary from %d file(s)]", scalar(@paths);
		    printf "\nTotal Read Count (R1 + R2)\t%s", $n;                                                                                # total sequencing read number
		    printf "\nQC-failed Read Count\t%s", $count{failed};                  # QC-failed reads in flagstat
		    printf "\nUnique Read Pairs\t%s\t%s (%%)", $count{"unique\tprimary\tfirst"} + $count{"unique\tunmapped\tfirst"}, ($count{"unique\tprimary\tfirst"} + $count{"unique\tunmapped\tfirst"}) / $n * 100 * 2;
		    printf "\nTotal Mapped Reads\t%s\t%s (%%)", $count{"duplicate\tprimary"} + $count{"unique\tprimary"}, ($count{"duplicate\tprimary"} + $count{"unique\tprimary"}) / $n * 100;
		    printf "\nNon-Mapped Reads\t%s\t%s (%%)", $count{"duplicate\tunmapped"} + $count{"unique\tunmapped"}, ($count{"duplicate\tunmapped"} + $count{"unique\tunmapped"}) / $n * 100;
		    printf "\nUnique Mapped Reads\t%s\t%s (%%)", $count{"unique\tprimary"}, $count{"unique\tprimary"} / $n * 100;
		    printf "\nMapped Read Duplication\t%s\t%s (%%)", $count{"duplicate\tprimary"}, $count{"duplicate\tprimary"} / $n * 100;
		    printf "\nStrand ratio (forward, reverse, reverse/forward of unique mapped)\t%s\t%s\t%s", $count{"unique\tprimary\tforward"}, $count{"unique\tprimary\treverse"}, $count{"unique\tprimary\treverse"} / $count{"unique\tprimary\tforward"};
		    
		    # Sequence length statistics
		    printf "\n\n\n[Sequencing read alignment statistics]\n";
		    printhash2(%count);
		    
		    
		    # the end of the main subroutine
		}



		## ------------ Library of subroutines --------------------------------------- #


		sub printhash2 {
		    # prints a hash
		    my (%hash) = @_;

		    # local variables
		    
		    foreach my $key (sort keys %hash)
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
        default: "$(inputs.input_files[0].basename).txt"
stdout: $(inputs.output_name)
outputs:
    aligned_stats:
        type: stdout

