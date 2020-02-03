#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Generate vaf report for concordance analysis"
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:bionic"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'eval_vaf_report.pl'
        entry: |
            #!/usr/bin/perl

            use strict;
            use warnings;

            use feature qw(say);

            die("wrong number of arguments, 10 are needed") unless scalar(@ARGV) == 10;
            my ($out_dir, $vcf, $b_n_snv_br, $b_n_indel_br, $b_t_snv_br, $b_t_indel_br, $q_n_snv_br, $q_n_indel_br, $q_t_snv_br, $q_t_indel_br) = @ARGV;

            my $b_n_snv_bamrc_info   = get_bamrc_info($b_n_snv_br);
            my $b_n_indel_bamrc_info = get_bamrc_info($b_n_indel_br);
            my $b_t_snv_bamrc_info   = get_bamrc_info($b_t_snv_br);
            my $b_t_indel_bamrc_info = get_bamrc_info($b_t_indel_br);

            my $q_n_snv_bamrc_info   = get_bamrc_info($q_n_snv_br);
            my $q_n_indel_bamrc_info = get_bamrc_info($q_n_indel_br);
            my $q_t_snv_bamrc_info   = get_bamrc_info($q_t_snv_br);
            my $q_t_indel_bamrc_info = get_bamrc_info($q_t_indel_br);
            
            my $vaf_out = $out_dir.'/vaf_report.out';
            open(my $out_fh, ">$vaf_out") or die "fail to write to $vaf_out";

            say $out_fh join "\t", 'Chr', 'Pos', 'Ref', 'Alt', 'Base', 'Query', 'Base_Normal_Dp', 'Base_Normal_Vaf', 'Base_Tumor_Dp', 'Base_Tumor_Vaf',  'Query_Normal_Dp', 'Query_Normal_Vaf', 'Query_Tumor_Dp', 'Query_Tumor_Vaf';

            open(my $fh, "/bin/gunzip -c $vcf |") or die("couldn't open $vcf to read");

            while (<$fh>) {
                next if $_ =~ /^#/;
                my @columns = split /\t/, $_;
                my ($chr, $pos, $ref, $alt, $info) = map{$columns[$_]}qw(0 1 3 4 7);
    
                my ($base, $query);
                if ($info =~ /source=Intersection/) {
                    ($base, $query) = (1, 1);
                }
                elsif ($info =~ /source=base/) {
                    ($base, $query) = (1, 0);
                }
                elsif ($info =~ /source=query/) {
                    ($base, $query) = (0, 1);
                }
                elsif ($info =~ /source=(FilteredInAll|filterInquery-base|query-filterInbase)/) {
                    ($base, $query) = ('-', '-');
                }
                else {
                    die "Line $_ does not have source info";
                }
    
                my @alts = split /,/, $alt;
                my ($b_n_AF, $b_n_dp) = get_AF_dp($chr, $pos, $ref, $alts[0], $b_n_snv_bamrc_info, $b_n_indel_bamrc_info);
                my ($b_t_AF, $b_t_dp) = get_AF_dp($chr, $pos, $ref, $alts[0], $b_t_snv_bamrc_info, $b_t_indel_bamrc_info);
                my ($q_n_AF, $q_n_dp) = get_AF_dp($chr, $pos, $ref, $alts[0], $q_n_snv_bamrc_info, $q_n_indel_bamrc_info);
                my ($q_t_AF, $q_t_dp) = get_AF_dp($chr, $pos, $ref, $alts[0], $q_t_snv_bamrc_info, $q_t_indel_bamrc_info);
                
                say $out_fh join "\t", $chr, $pos, $ref, $alt, $base, $query, $b_n_dp, $b_n_AF, $b_t_dp, $b_t_AF, $q_n_dp, $q_n_AF, $q_t_dp, $q_t_AF;
            }
            $fh->close;

            sub get_bamrc_info {
                my $bamrc = shift;
                my $bamrc_info;
    
                open(my $bamrc_fh, $bamrc) or die "failed to open $bamrc for read";
                while (<$bamrc_fh>) {
                    chomp;
                    my ($chr, $pos, $ref, $depth, $rest) = $_ =~ /^(.+)\t(\d+)\t([AaCcGgTtRrYyKkMmSsWwBbDdHhVvNn])\t(\d+)\t(.+)$/x; 
                    my @allele_info = split /\t/, $rest;
                    my $allele_ct;
                    for my $allele_str (@allele_info) {
                        my @info = split /:/, $allele_str;
                        $allele_ct->{$info[0]} = $info[1];
                    }
                    $bamrc_info->{$chr.'__'.$pos} = {
                        ref => $ref,
                        depth => $depth,
                        allele_ct => $allele_ct,
                    };
                }
                close $bamrc_fh;
                return $bamrc_info;
            }

            sub get_AF_dp {
                my ($chr, $pos, $ref, $alt, $snv_br_info, $indel_br_info) = @_;
                my ($AF, $dp) = ('NA', 'NA');
                
                if (length($ref) == 1 and length($alt) == 1) { #SNV
                    my $id = $chr .'__'.$pos;
                    my $br_info = $snv_br_info->{$id};
                    if ($br_info) {
                        $AF = sprintf("%.5f", $br_info->{allele_ct}->{$alt}/$br_info->{depth});
                        $dp = $br_info->{depth};
                    }
                }
                else {
                    my ($indel_pos, $bamrc_str);
                    if (length($ref) > length($alt)) { #DEL
                        $indel_pos = $pos + 1;
                        $bamrc_str = '-'.substr($ref, length($alt));
                    }
                    else { #INS
                        $indel_pos = $pos;
                        $bamrc_str = '+'.substr($alt, length($ref));
                    }
                    my $id = $chr.'__'.$indel_pos;
                    my $br_info = $indel_br_info->{$id};
                    if ($br_info) {
                        my $alt_ct = exists $br_info->{allele_ct}->{$bamrc_str} ? $br_info->{allele_ct}->{$bamrc_str} : 0;
                        $AF = sprintf("%.5f", $alt_ct/$br_info->{depth});
                        $dp = $br_info->{depth};
                    }
                }
                return ($AF, $dp);
            }

arguments: [
    "/usr/bin/perl", "eval_vaf_report.pl", $(runtime.outdir)
]
inputs:
    combined_vcf:
        type: File
        inputBinding:
            position: 1
    base_normal_snv_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 2
    base_normal_indel_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 3
    base_tumor_snv_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 4
    base_tumor_indel_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 5
    query_normal_snv_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 6
    query_normal_indel_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 7
    query_tumor_snv_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 8
    query_tumor_indel_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 9
outputs:
    out_file:
        type: File
        outputBinding:
            glob: "vaf_report.out"

