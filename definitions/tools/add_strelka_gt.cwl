#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add GT tags"
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:bionic"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'add_strelka_gt.pl'
        entry: |
            #!/usr/bin/perl

            use strict;
            use warnings;

            use feature qw(say);

            die("wrong number of arguments") unless scalar(@ARGV) == 2;
            my ($strelka_vcf, $outdir) = @ARGV;

            open(my $strelka_vcf_fh, '-|', '/bin/gunzip', '-c', $strelka_vcf) 
                or die("couldn't open $strelka_vcf to read");
            open(my $add_gt_fh, ">", "$outdir/add_gt.vcf") 
                or die("couldn't open add_gt.vcf for write");

            while (<$strelka_vcf_fh>) {
                chomp;
                if (/^##/) {
                    say $add_gt_fh $_;
                }
                elsif (/^#/) { #COLUMN HEADER
                    say $add_gt_fh '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
                    say $add_gt_fh $_;
                }
                else {
                    my @columns = split /\t/, $_;
                    my ($ref, $alt, $info) = map{$columns[$_]}(3, 4, 7);
                    my @alts = split /,/, $alt;

                    my ($n_gt, $t_gt);

                    if (length($ref) == 1 and length($alts[0]) == 1) {
                        my ($n_gt_info, $n_gt_str, $t_gt_str) = $info =~ /NT=(\S+?);QSS.*SGT=(\S+?)\->(\S+?);/;
                        unshift @alts, $ref;

                        my %ids;
                        my $id = 0;

                        for my $base (@alts) {
                            $ids{$base} = $id;
                            $id++;
                        }

                        $n_gt = $n_gt_info eq 'ref' ? '0/0' : parse_gt($n_gt_str, \%ids);
                        $t_gt = parse_gt($t_gt_str, \%ids);
                    }
                    else {#INDEL
                        my ($n_gt_info, $t_gt_info) = $info =~ /;NT=(\S+?);.*SGT.*\->(\S+?);/;

                        my %gt_info = (
                            ref => '0/0',
                            het => '0/1',
                            hom => '1/1',
                            conflict => './.',
                        );

                        $n_gt = $gt_info{$n_gt_info};
                        $t_gt = $gt_info{$t_gt_info};
                    }

                    $columns[8]  = 'GT:'.$columns[8];
                    $columns[9]  = $n_gt . ':' . $columns[9];
                    $columns[10] = $t_gt . ':' . $columns[10];

                    my $new_str = join "\t", @columns;
                    say $add_gt_fh $new_str;
                }
            }

            close($strelka_vcf_fh);
            close($add_gt_fh);


            sub parse_gt {
                my ($gt_str, $ids) = @_;
                my @gt_ids = map{$ids->{$_}}(split //, $gt_str);
                return join '/', sort @gt_ids;
            }

            #SNV example
            #1       10231   .       C       A       .       QSS_ref NT=ref;QSS=1;QSS_NT=1;SGT=AC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    32:4:8:0:0,3:28,60:0,0:0,1      84:6:69:0:7,21:71,192:0,0:0,1
            #INDEL example
            ##1     965051  .       ATGTGTG A       .       QSI_ref IC=5;IHP=2;NT=ref;QSI=1;QSI_NT=1;RC=8;RU=TG;SGT=het->het;SOMATIC;TQSI=1;TQSI_NT=1       DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50   8:8:6,6:0,0:2,4:10.3:0.00:0.00  18:18:8,8:5,6:5,8:21:0.25:0.00

arguments: [
    "/usr/bin/perl", "add_strelka_gt.pl",
    $(inputs.vcf.path), $(runtime.outdir)
]
inputs:
    vcf:
        type: File
outputs:
    processed_vcf:
        type: File
        outputBinding:
            glob: "add_gt.vcf"

