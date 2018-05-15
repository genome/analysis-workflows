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
for my $file_input ('omni_vcf') {
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

for my $other_input ('picard_metric_accumulation_level') {
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

my $target_region_set_name;
my $multiple_trsn = 0;

for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    push @bams, {class => 'File', path => $id->bam_path};

    unless ($target_region_set_name) {
        $target_region_set_name = $id->target_region_set_name;
    } else {
        if ($id->target_region_set_name ne $target_region_set_name) {
            $multiple_trsn = 1;
        }
    }
}

$inputs{bam} = \@bams;

my ($roi, @extra_roi) = grep { $_->name eq 'region_of_interest_set_name' } @inputs;
if (@extra_roi) {
    $build->fatal_message('multiple inputs found for reference');
}
my $fl;
if ($roi) {
    $fl = Genome::FeatureList->get(name => $roi)
        or $build->fatal_message('no ROI found for input');
} elsif ($multiple_trsn) {
    $build->fatal_message('Multiple target_region_set_name found on instrument data. Please specify a "region_of_interest_set_name" input explicitly.');
} elsif ($target_region_set_name) {
    $fl = Genome::FeatureList->get(name => $target_region_set_name)
        or $build->fatal_message('no ROI found for instrument data target region set');
} else {
    $build->fatal_message('No target_region_set_name found on instrument data.  Please specify a "region_of_interest_set_name" input explicitly.');
}

my $sr_users = Genome::SoftwareResult::User->user_hash_for_build($build);

my @intervals_params = (
    feature_list => $fl,
    reference_build => $ref_build,
    users => $sr_users,
    merge => 1,
);

my $target_interval_sr = Genome::FeatureList::IntervalList->get_with_lock(
    @intervals_params,
    track_name => 'target_region',
);
my $bait_interval_sr = Genome::FeatureList::IntervalList->get_with_lock(
    @intervals_params,
    track_name => 'tiled_region',
);
unless($target_interval_sr and $bait_interval_sr) {
    $build->fatal_message('No interval list results found.  Please run `genome feature-list dump-interval-list --feature-list %s --reference %s --merge` twice, once with each of `--track "target_region"` and `--track "tiled_region"` options.', $fl->id, $ref_build->id);
}

$inputs{target_intervals} = {
        class => 'File',
        path => $target_interval_sr->interval_list,
};

$inputs{bait_intervals} = {
        class => 'File',
        path => $bait_interval_sr->interval_list,
};
my @read_structure_inputs = grep { $_->name eq 'read_structure' } @inputs;
my @read_structure = map $_->value_id, @read_structure_inputs;

unless(@read_structure) {
    push @read_structure, ('8M+T') x 2;
}

$inputs{read_structure} = \@read_structure;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, \%inputs);

