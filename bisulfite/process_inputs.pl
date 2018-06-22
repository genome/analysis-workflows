#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Spec;
use YAML::XS;

my $build_id = $ARGV[0]
    or die 'no build id';

my $build = Genome::Model::Build->get($build_id)
    or die 'no build for id';

my @inputs = $build->inputs;


my %inputs;
for my $file_input ('reference_sizes', 'trimming_adapters') {
    my ($found, @extra) = grep { $_->name eq $file_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $file_input;
    }

    $inputs{$file_input} = {
        class => 'File',
        path => $found->value_id,
    };
}

for my $string_input ('reference_index', 'trimming_adapter_trim_end','trimming_adapter_min_overlap','trimming_max_uncalled','trimming_min_readlength') {
    my ($found, @extra) = grep { $_->name eq $string_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $string_input;
    }

    $inputs{$string_input} = $found->value_id;
}
#$inputs{reference} = $inputs{'reference_index'};

my (@bams, @rg_ids, @rg_fields);

my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    push @bams, {class => 'File', path => $id->bam_path};
    push @rg_ids, $id->id;

    my $pu = join('.', $id->flow_cell_id, $id->lane, ( $id->can('index_sequence')? $id->index_sequence : () ));
    my $sm = $id->sample->name;
    my $lb = $id->library->name;
    my $pl = 'Illumina';
    my $cn = 'WUGSC';
    my $rgid = $id->id;

    push @rg_fields, join("\t", '@RG', "ID:$rgid", "PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn");
}

$inputs{instrument_data_bams} = \@bams;
$inputs{read_group_id} = \@rg_ids;
$inputs{sample_name} = $build->subject->name;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, \%inputs);
