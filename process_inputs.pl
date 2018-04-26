#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Array::Compare;
use File::Spec;
use YAML::XS;

my $build_id = $ARGV[0]
    or die 'no build id';

my $build = Genome::Model::Build->get($build_id)
    or die 'no build for id';

my @inputs = $build->inputs;

my %inputs;
for my $file_input ('dbsnp', 'known_indels', 'mills', 'bait_intervals', 'target_intervals', 'omni_vcf', 'docm_vcf', 'synonyms_file', 'custom_gnomad_vcf') {
    my ($found, @extra) = grep { $_->name eq $file_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $file_input;
    }

    $inputs{$file_input} = {
        class => 'File',
        path => $found->value_id,
    };
}

$inputs{sample_name} = $build->subject->name;

for my $numeric_input ('varscan_strand_filter', 'varscan_min_coverage', 'varscan_min_var_freq', 'varscan_p_value', 'varscan_min_reads') {
    my ($found, @extra) = grep { $_->name eq $numeric_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $numeric_input;
    }

    $inputs{$numeric_input} = $found->value_id + 0;
}
for my $other_input ('picard_metric_accumulation_level', 'vep_cache_dir') {
    my ($found, @extra) = grep { $_->name eq $other_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $other_input;
    }

    $inputs{$other_input} = $found->value_id;
}

my ($ref, @extra) = grep { $_->name eq 'reference_build' } @inputs;
if (@extra) {
    die 'multiple inputs found for reference';
}

my $ref_build = $ref->value_class_name->get($ref->value_id)
    or die 'no reference found for input';
$inputs{reference} = $ref_build->full_consensus_path('fa');


my (@bams, @rg_fields);

my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
my @ids = map $_->value_id, @instrument_data_inputs;

for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    push @bams, {class => 'File', path => $id->bam_path};

    my $pu = join('.', $id->flow_cell_id, $id->lane, ( $id->can('index_sequence')? $id->index_sequence : () ));
    my $sm = $id->sample->name;
    my $lb = $id->library->name;
    my $pl = 'Illumina';
    my $cn = 'WUGSC';
    my $rgid = $id->id;

    push @rg_fields, join("\t", '@RG', "ID:$rgid", "PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn");
}

$inputs{bams} = \@bams;
$inputs{readgroups} = \@rg_fields;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, \%inputs);

