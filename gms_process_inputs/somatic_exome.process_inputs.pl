#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use YAML::XS;

my $build_id = $ARGV[0]
    or die 'no build id';

my $build = Genome::Model::Build->get($build_id)
    or die 'no build for id';

{
    package InputProcessor;
    class InputProcessor {
        is => 'Genome::Model::Build::CwlPipeline::InputProcessor',
        has_simple_input => [
            picard_metric_accumulation_level => { input_type => 'Text' },
            vep_cache_dir => { input_type => 'Text' },
            mills => { input_type => 'File' },
            known_indels => { input_type => 'File' },
            dbsnp_vcf => { input_type => 'File' },
            docm_vcf => { input_type => 'File' },
            synonyms_file => { input_type => 'File' },
            omni_vcf => { input_type => 'File' },
            cosmic_vcf => { input_type => 'File' },
            panel_of_normals_vcf => { input_type => 'File' },
            custom_gnomad_vcf => { input_type => 'File' },
            custom_clinvar_vcf => { input_type => 'File' },
            vep_assembly => { input_type => 'Text' },
            manta_call_regions => { input_type => 'File' },
            manta_non_wgs => { input_type => 'Text' },
            manta_output_contigs => { input_type => 'Text' },
            somalier_vcf => { input_type => 'File' },
        ],
    };
}
my @inputs = $build->inputs;
my $input_processor = InputProcessor->get($build->id);
my $inputs = $input_processor->simple_inputs;
$inputs->{mutect_scatter_count} = 50;
$inputs->{bqsr_intervals} = [map 'chr'.$_, 1..22];
my ($ref, @extra) = grep { $_->name eq 'reference_build' } @inputs;
if (@extra) {
    die 'multiple inputs found for reference';
}

my $ref_build = $ref->value_class_name->get($ref->value_id)
    or die 'no reference found for input';
$inputs->{reference} = $ref_build->full_consensus_path('fa');
my (@normal_bams, @normal_rg_fields, @tumor_bams, @tumor_rg_fields);

my ($tumor_input, @extra) = grep { $_->name eq 'tumor_sample' } @inputs;
if (@extra) {
    die 'multiple inputs found for tumor sample';
}


my $tumor_sample = $tumor_input->value_class_name->get($tumor_input->value_id)
    or die 'no tumor found for input';

my ($normal_input, @extra) = grep { $_->name eq 'normal_sample' } @inputs;
if (@extra) {
    die 'multiple inputs found for normal sample';
}

my $normal_sample = $normal_input->value_class_name->get($normal_input->value_id)
    or die 'no normal found for input';

my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
my @ids = map $_->value_id, @instrument_data_inputs;

my $target_region_set_name;
my $multiple_trsn = 0;

for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';
    
    my $pu = join('.', $id->flow_cell_id, $id->lane, ( $id->can('index_sequence')? $id->index_sequence : () ));
    my $sm = $id->sample->name;
    my $lb = $id->library->name;
    my $pl = 'Illumina';
    my $cn = 'WUGSC';
    my $rgid = $id->id;

    my $bam_path = $id->bam_path
      or die 'no bam path found for instrument data ' . $id->id;
    
    if($id->sample eq $tumor_sample){
      push @tumor_bams, {class => 'File', path => $bam_path};
      push @tumor_rg_fields, join("\t", '@RG', "ID:$rgid", "PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn");
    } elsif($id->sample eq $normal_sample) { 
      push @normal_bams, {class => 'File', path => $bam_path};
      push @normal_rg_fields, join("\t", '@RG', "ID:$rgid", "PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn");
    } else {
       die "no sample match found:" . $id->id;
    }

    unless ($target_region_set_name) {
        $target_region_set_name = $id->target_region_set_name;
    } else {
        if ($id->target_region_set_name ne $target_region_set_name) {
            $multiple_trsn = 1;
        }
    }
}

$inputs->{tumor_bams} = \@tumor_bams;
$inputs->{tumor_readgroups} = \@tumor_rg_fields;

$inputs->{normal_bams} = \@normal_bams;
$inputs->{normal_readgroups} = \@normal_rg_fields;

my ($roi_input, @extra_roi) = grep { $_->name eq 'region_of_interest_set_name' } @inputs;
if (@extra_roi) {
    $build->fatal_message('multiple inputs found for reference');
}
my $fl;

if ($roi_input) {
    $fl = Genome::FeatureList->get(name => $roi_input->value_id)
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

$inputs->{target_intervals} = {
        class => 'File',
        path => $target_interval_sr->interval_list,
};

$inputs->{bait_intervals} = {
        class => 'File',
        path => $bait_interval_sr->interval_list,
};

$inputs->{interval_list} ={
        class => 'File',
        path => $target_interval_sr->interval_list,
};

$inputs->{per_base_intervals} = [
    { label => 'clinvar', file => '/gscmnt/gc2560/core/model_data/interval-list/01f4fae3699646c3af2fa47853da7a8c/06a82ecf9c434b7ab03d82e59eaa28c8.interval_list' },
];
$inputs->{per_target_intervals} = [
    { label => 'acmg_genes', file => '/gscmnt/gc2560/core/model_data/interval-list/db8c25932fd94d2a8a073a2e20449878/a35b64d628b94df194040032d53b5616.interval_list' },
];
$inputs->{summary_intervals} = [];

$inputs->{variants_to_table_fields} = [qw(CHROM POS REF ALT set)];
$inputs->{variants_to_table_genotype_fields} = [qw(GT AD AF DP)];
$inputs->{vep_to_table_fields} = [qw(Consequence SYMBOL Feature_type Feature HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons HGNC_ID Existing_variation gnomADe_AF CLIN_SIG SOMATIC PHENO clinvar_CLINSIGN clinvar_PHENOTYPE clinvar_SCORE clinvar_RCVACC clinvar_TESTEDINGTR clinvar_PHENOTYPELIST clinvar_NUMSUBMIT clinvar_GUIDELINES)];

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, $inputs);
