#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bam-readcount"

baseCommand: ["/usr/bin/perl", "helper.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/bam-readcount-cwl:0.7.4"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'helper.pl'
        entry: |
            use warnings;
            use strict;

            use IO::File;
            use File::Temp;

            my $file = $ARGV[0];
            my $sample = $ARGV[1];
            my $bam_file = $ARGV[2];
            my $ref_fasta = $ARGV[3];
            my $snp_output = $ARGV[4];
            my $indel_output = $ARGV[5];

            my $input = IO::File->new($file);
            my @vcf_header = parse_vcf_header($input);
            my $header_line = $vcf_header[-1];
            chomp $header_line;
            my @header_fields = split \"\\t\", $header_line;

            my @vcf_lines; #20
            my %rc_for_indel;
            my %rc_for_snp;
            while (my $entry = $input->getline) {
                push @vcf_lines, $entry;
                chomp $entry;

                my %fields;
                @fields{@header_fields} = split \"\\t\", $entry;

                my $filter_sample = $fields{$sample};
                unless ($filter_sample) {
                    die "Unable to find field for $sample";
                }

                my @sample_fields = split ':', $filter_sample;
                unless (@sample_fields) {
                    die "Unable to parse field for $sample";
                }

                my $index = 0; #40
                my %format_keys = map { $_ => $sample_fields[$index++] } split ':', $fields{FORMAT};
                #these are in order ACGT
                my @alleles = ($fields{REF}, split ',', $fields{ALT});
                my %gt_alleles = map {$_ => 1} grep { $_ > 0 } grep { $_ ne '.' } split '/', $format_keys{GT};
                my @used_alleles;
                for my $allele_index (keys %gt_alleles) {
                    push @used_alleles, $alleles[$allele_index];
                }
                my ($var) = sort @used_alleles;
                $var = q{} unless defined $var;
                $var = uc($var);
                my $ref = uc($fields{REF});
                my $chrom = $fields{'#CHROM'};
                my $pos = $fields{POS};

                if (length($ref) > 1 || length($var) > 1) {
                    #it's an indel or mnp
                    if (length($ref) == length($var)) {
                        die "MNPs unsupported";
                    } #60
                    elsif (length($ref) > length($var)) {
                        #it's a deletion
                        $pos += 1;
                        ($ref, $var) = ($var, $ref);
                        $ref = substr($var, 1, 1);
                        $var = "-" . substr($var, 1);
                    }
                    else {
                        #it's an insertion
                        substr($var, 0, 1, "+");
                    }
                    $rc_for_indel{$chrom}{$pos}{$ref}{$var} = $entry;
                }
                else {
                    #it's a snp
                    $rc_for_snp{$chrom}{$pos}{$ref}{$var} = $entry;
                }
            }

            if (%rc_for_indel) { #80
                write_bam_readcount_file(\\%rc_for_indel, $bam_file, $ref_fasta, $indel_output, '-i')
            }
            if (%rc_for_snp) {
                write_bam_readcount_file(\\%rc_for_snp, $bam_file, $ref_fasta, $snp_output)
            }

            sub write_bam_readcount_file {
                my ($hash, $bam_file, $ref_fasta, $output_file, $optional_param) = @_;
                $optional_param ||= '';
                my ($list_fh, $list_name) = File::Temp::tempfile();
                generate_region_list($hash, $list_fh);
                $list_fh->close();
                `/usr/bin/bam-readcount -f $ref_fasta -l $list_name -w 0 -b 20 $optional_param $bam_file > $output_file`
            }

            sub generate_region_list {
                my ($hash, $region_fh) = @_; #input_fh should be a filehandle to the VCF
                for my $chr (keys %$hash) {
                    for my $pos (sort { $a <=> $b } keys %{$hash->{$chr}}) { #100
                        print $region_fh join("\t", $chr, $pos, $pos) . "\n";
                    }
                }
            }

            sub parse_vcf_header {
                my $input_fh = shift;
                my @header;
                my $header_end = 0;
                while (!$header_end) {
                    my $line = $input_fh->getline;
                    if ($line =~ m/^##/) {
                        push @header, $line;
                    }
                    elsif ($line =~ m/^#/) {
                        push @header, $line;
                        $header_end = 1
                    }
                    else {
                        die "Malformed header:\n" . join("\n", @header)
                    }
                }
                return @header;
            }

arguments:
    - position: 5
      valueFrom: $(runtime.outdir)/snp.bam_readcount
    - position: 6
      valueFrom: $(runtime.outdir)/indel.bam_readcount
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    sample:
        type: string
        inputBinding:
            position: 2
    bam:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [^.bai]
    reference_fasta:
        type: File
        inputBinding:
            position: 4
        secondaryFiles: [.fai]
outputs:
    snp_bam_readcount:
        type: File?
        outputBinding:
            glob: snp.bam_readcount
    indel_bam_readcount:
        type: File?
        outputBinding:
            glob: indel.bam_readcount
