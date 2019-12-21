#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Make CLE gold evaluation output and vaf report"
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:bionic"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'eval_cle_gold.pl'
        entry: |
            #!/usr/bin/perl

            use strict;
            use warnings;

            use feature qw(say);

            die("wrong number of arguments, 10 are needed") unless scalar(@ARGV) == 10;
            my ($out_dir, $sompy_out, $tn_bed, $tn_query_snv_bed, $tn_query_indel_bed, $vcf, $n_snv_br, $n_indel_br, $t_snv_br, $t_indel_br) = @ARGV;

            my $tn_size = bed_size($tn_bed);
            my $tn_query_snv_size = bed_size($tn_query_snv_bed);
            my $tn_query_indel_size = bed_size($tn_query_indel_bed);

            my $snv_specificity   = sprintf("%.9f", ($tn_size-$tn_query_snv_size)/$tn_size);
            my $indel_specificity = sprintf("%.9f", ($tn_size-$tn_query_indel_size)/$tn_size);
            open(my $sompy_fh, $sompy_out) or die("couldn't open $sompy_out to read");
            open(my $eval_out_fh, ">$out_dir/gold_eval.out") or die("couldn't open gold_eval.out for write");
          
            while (<$sompy_fh>) {
               next if /^\s+$|\s+records\s+/;
               chomp;
               if (/^\s+type\s+/) {
                   say $eval_out_fh $_."  specificity";
               }
               elsif (/\s+indels\s+/) {
                   say $eval_out_fh $_."  ".$indel_specificity;
               }
               elsif (/\s+SNVs\s+/) {
                   say $eval_out_fh $_."  ".$snv_specificity;
               }
            }
            close $sompy_fh;
            close $eval_out_fh;

            my $t_snv_bamrc_info   = get_bamrc_info($t_snv_br);
            my $t_indel_bamrc_info = get_bamrc_info($t_indel_br);
            my $n_snv_bamrc_info   = get_bamrc_info($n_snv_br);
            my $n_indel_bamrc_info = get_bamrc_info($n_indel_br);

            my $vaf_out = $out_dir. '/vaf_report.out';
            open(my $out_fh, ">$vaf_out") or die "fail to write to $vaf_out";

            say $out_fh join "\t", 'Chr', 'Pos', 'Ref', 'Alt', 'Gold', 'Query', 'Normal_Dp', 'Normal_Vaf', 'Tumor_Dp', 'Tumor_Vaf';

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
                my ($n_AF, $n_dp) = get_AF_dp($chr, $pos, $ref, $alts[0], $n_snv_bamrc_info, $n_indel_bamrc_info);
                my ($t_AF, $t_dp) = get_AF_dp($chr, $pos, $ref, $alts[0], $t_snv_bamrc_info, $t_indel_bamrc_info);
                say $out_fh join "\t", $chr, $pos, $ref, $alt, $base, $query, $n_dp, $n_AF, $t_dp, $t_AF;
            }
            $fh->close;

            sub bed_size {
                my $bed = shift;
                my $size = 0;
                open(my $fh, $bed) or die "could not open $bed for read";
                while (<$fh>) {
                    my (undef, $start, $end) = split /\t/, $_;
                    $size += $end - $start;
                }
                close $fh;
                return $size;
            }

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
    "/usr/bin/perl", "eval_cle_gold.pl", $(runtime.outdir)
]
inputs:
    sompy_out:
        type: File
        inputBinding:
            position: 1
    true_negative_bed:
        type: File
        inputBinding:
            position: 2
    tn_x_query_snv:
        type: File
        inputBinding:
            position: 3
    tn_x_query_indel:
        type: File
        inputBinding:
            position: 4
    combined_vcf:
        type: File
        inputBinding:
            position: 5
    normal_snv_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 6
    normal_indel_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 7
    tumor_snv_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 8
    tumor_indel_bam_readcount_tsv:
        type: File
        inputBinding:
            position: 9
outputs:
    eval_out:
        type: File
        outputBinding:
            glob: "gold_eval.out"
    vaf_report:
        type: File
        outputBinding:
            glob: "vaf_report.out"

