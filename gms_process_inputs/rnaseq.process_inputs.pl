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
for my $file_input ('reference_index', 'trimming_adapters', 'reference_annotation', 'kallisto_index', 'gene_transcript_lookup_table','refFlat','ribosomal_intervals') {
    my ($found, @extra) = grep { $_->name eq $file_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $file_input;
    }

    $inputs{$file_input} = {
        class => 'File',
        path => $found->value_id,
    };
}

for my $number_input ('trimming_adapter_min_overlap', 'trimming_max_uncalled', 'trimming_min_readlength') {
    my ($found, @extra) = grep { $_->name eq $number_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $number_input;
    }

    $inputs{$number_input} = int($found->value_id);
}

for my $value_input ('trimming_adapter_trim_end','strand') {
    my ($found, @extra) = grep { $_->name eq $value_input } @inputs;

    if (@extra) {
        die 'multiple inputs found for ' . $value_input;
    }

    $inputs{$value_input} = $found->value_id;
}


my (@bams, @rg_ids, @rg_fields);

my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    my $bam_path = $id->bam_path
        or die 'no bam_path found for instrument data ' . $id->id;

    push @bams, {class => 'File', path => $id->bam_path};
    push @rg_ids, $id->id;

    my $pu = join('.', $id->flow_cell_id, $id->lane, ( $id->can('index_sequence')? $id->index_sequence : () ));
    my $sm = $id->sample->name;
    my $lb = $id->library->name;
    my $pl = 'Illumina';
    my $cn = 'WUGSC';

    push @rg_fields, ["PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn"];
}

$inputs{instrument_data_bams} = \@bams;
$inputs{read_group_id} = \@rg_ids;
$inputs{read_group_fields} = \@rg_fields;
$inputs{sample_name} = $build->subject->name;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, \%inputs);
