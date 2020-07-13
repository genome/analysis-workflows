#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use YAML::XS;
use Genome::Utility::Text qw( sanitize_string_for_filesystem );

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
            synonyms_file => { input_type => 'File' },
            omni_vcf => { input_type => 'File' },
            idt_bed => { input_type => 'File' },
            intervals_dir => { input_type => 'File' },
            bqsr_interval_list => { input_type => 'Text' },
            vep_ensembl_assembly => { input_type => 'Text' },
            vep_ensembl_version => { input_type => 'Text' },
            vep_ensembl_species => { input_type => 'Text' },
        ],
    };
}

my $intervals_file = '/gscmnt/gc2698/jin810/references/split_regions/split_regions.txt';
my @lines = Genome::Sys->read_file($intervals_file); 
chomp(@lines);
#my ($intervals_input, @extra) = $build->inputs(name => 'intervals_file'); 
#$build->fatal_message('too many inputs for intervals_file') if @extra;
#$build->fatal_message('missing intervals_file input') unless ($intervals_input);  #if needs to be required
#my $intervals_file = $intervals_input->value_id;

my @inputs = $build->inputs;
my $input_processor = InputProcessor->get($build->id);
my $inputs = $input_processor->simple_inputs;
$inputs->{bqsr_intervals} = [ map { 'chr'.$_ } (1..22, 'X', 'Y') ];
$inputs->{intervals} = [map { [$_] } @lines];
my ($ref, @extra) = grep { $_->name eq 'reference_build' } @inputs;
if (@extra) {
    die 'multiple inputs found for reference';
}
my $sr_users = Genome::SoftwareResult::User->user_hash_for_build($build);

my $ref_build = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(users => $sr_users, reference_build_id => $ref->value_id, aligner_name => 'bwamem', aligner_version => '0.7.15')
    or die 'no bwamem 0.7.15 aligner index found for reference input';

$inputs->{reference} = $ref_build->full_consensus_path('fa');



$inputs->{sample_name} = $build->subject->name;


$inputs->{emit_reference_confidence} = 'BP_RESOLUTION';
$inputs->{gvcf_gq_bands} = [];
#$inputs->{intervals} = [ map { ['chr'.$_] } (1..22, 'X', 'Y') ];

my @sequence;
my $target_region_set_name;
my $multiple_trsn = 0;

my @instrument_data_inputs = grep { $_->name eq 'instrument_data' } @inputs;
my @ids = map $_->value_id, @instrument_data_inputs;

for my $input (@instrument_data_inputs) {
    my $id = $input->value_class_name->get($input->value_id)
        or die 'no instrument data found for input';

    my $pu = join('.', $id->flow_cell_id, $id->lane, ( $id->can('index_sequence')? $id->index_sequence : () ));
    my $sm = $id->sample->name;
    my $lb = $id->library->name;
    my $pl = 'Illumina';
    my $cn = 'WUGSC';
    my $rgid = $id->id;

    my $sequence = {};
    $sequence->{readgroup} = join("\t", '@RG', "ID:$rgid", "PU:$pu", "SM:$sm", "LB:$lb", "PL:$pl", "CN:$cn");

    if (my $bam_path = $id->bam_path) {
        $sequence->{sequence}{bam} = {class => 'File', path => $bam_path};
    } elsif (my @fastqs = sort glob(File::Spec->join($id->disk_allocation->absolute_path, '*.fastq.gz'))) {
        unless (@fastqs == 2) {
            die "expected two fastqs but got " . scalar(@fastqs);
        }

        $sequence->{sequence}{fastq1} = { class => 'File', path => $fastqs[0] };
        $sequence->{sequence}{fastq2} = { class => 'File', path => $fastqs[1] };
    } else {
        die 'No FASTQs or BAM found for instrument data ' . $id->id;
    }
    push @sequence, $sequence;

    unless ($target_region_set_name) {
        $target_region_set_name = $id->target_region_set_name;
    } else {
        if ($id->target_region_set_name ne $target_region_set_name) {
            $multiple_trsn = 1;
        }
    }
}
$inputs->{sequence} = \@sequence;

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

my (@per_target_intervals, @per_base_intervals, @summary_intervals);

for my $type ('per_target', 'per_base') {
    my ($per_target_roi, @extra_pt_roi) = grep { $_->name eq "${type}_region_of_interest_set_name" } @inputs;
    if (@extra_pt_roi) {
        $build->fatal_message('multiple %s inputs found for reference', $type);
    }
    my $per_target_fl;
    if ($per_target_roi) {
        $per_target_fl = Genome::FeatureList->get(name => $per_target_roi->value_id)
            or $build->fatal_message('no %s ROI found for input', $type);
    } else {
        $build->fatal_message('Missing %s_region_of_interest_set_name', $type);
    }

    my @per_target_intervals_params = (
        feature_list => $per_target_fl,
        reference_build => $ref_build,
        users => $sr_users,
        merge => 1,
    );

    my $per_target_interval_sr = Genome::FeatureList::IntervalList->get_with_lock(
        @per_target_intervals_params,
        track_name => 'target_region',
    );
    unless($per_target_interval_sr) {
        $build->fatal_message('No interval list result found.  Please run `genome feature-list dump-interval-list --feature-list %s --reference %s --merge --track "target_region"`.', $per_target_fl->id, $ref_build->id);
    }

   my $sanitized_name = sanitize_string_for_filesystem($per_target_fl->name);
   push @per_target_intervals, { label => $sanitized_name, file => $per_target_interval_sr->interval_list };
   push @per_base_intervals, { label => $sanitized_name, file => $per_target_interval_sr->interval_list };
   push @summary_intervals, { label => $sanitized_name, file => $per_target_interval_sr->interval_list };
}

$inputs->{per_target_intervals} = \@per_target_intervals;
$inputs->{per_base_intervals} = \@per_base_intervals;
$inputs->{summary_intervals} = \@summary_intervals;

$inputs->{vep_plugins} = ['ExACpLI', 'LoFtool', 'SpliceRegion', 'dbNSFP,/gscmnt/gc2560/core/cwl/inputs/dbNSFP_scoring/dbNSFP4.0b2a.b38.gz,CADD_raw,CADD_raw_rankscore,CADD_phred,REVEL_score,REVEL_rankscore'];
$inputs->{variants_to_table_genotype_fields} = ['GT', 'AD', 'DP'];
$inputs->{variants_to_table_fields} = ['CHROM', 'POS', 'ID', 'REF', 'ALT'];
$inputs->{vep_to_table_fields} = ['Consequence', 'SYMBOL', 'Feature', 'HGVSc', 'HGVSp', 'ExACpLI', 'LoFtool', 'SpliceRegion', 'CADD_raw', 'REVEL_score', 'gnomADe_AF', 'Consequence', 'SYMBOL', 'Feature_type', 'Feature', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'HGNC_ID', 'Existing_variation', 'CLIN_SIG', 'PHENO', 'clinvar_CLINSIGN', 'clinvar_PHENOTYPE', 'clinvar_SCORE', 'clinvar_RCVACC', 'clinvar_TESTEDINGTR', 'clinvar_PHENOTYPELIST', 'clinvar_NUMSUBMIT', 'clinvar_GUIDELINES'];
#['Consequence', 'SYMBOL', 'Feature', 'HGVSc', 'HGVSp', 'ExACpLI', 'LoFtool', 'CADD_raw', 'REVEL_score', 'gnomADe_AF', 'Consequence', 'SYMBOL', 'Feature_type', 'Feature', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'HGNC_ID', 'Existing_variation', 'CLIN_SIG', 'PHENO'];

my @vep_custom_annotations;

my ($gnomad_file, @extra) = grep { $_->name eq 'custom_gnomad_vcf' } @inputs;
if (@extra) {
    die 'multiple inputs found for custom_gnomad_vcf';
}
if($gnomad_file) {
    my $vcf_path = $gnomad_file->value_id;
    my $custom_annotation->{method} = 'exact';
    $custom_annotation->{force_report_coordinates} = 'true';

    my $annotation_info->{data_format} = 'vcf';
    $annotation_info->{file} = { class => 'File', path => $vcf_path, secondaryFiles => [{class => 'File', path => $vcf_path . '.tbi' }]};
    $annotation_info->{name} = 'gnomadW'; # if changed should update the `vep_to_table_fields` input below
    $annotation_info->{gnomad_filter} = 'true';
    $annotation_info->{check_existing} = 'true';
    $annotation_info->{vcf_fields} = ['AF','AF_afr','AF_amr','AF_asj','AF_eas','AF_fin','AF_nfe','AF_oth','AF_sas'];
    $custom_annotation->{annotation} = $annotation_info;

    push @vep_custom_annotations, $custom_annotation;
}

my ($clinvar_file, @extra) = grep { $_->name eq 'custom_clinvar_vcf' } @inputs;
if (@extra) {
    die 'multiple inputs found for custom_clinvar_vcf';
}
if($clinvar_file) {
    my $vcf_path = $clinvar_file->value_id;
    my $custom_annotation->{method} = 'exact';
    $custom_annotation->{force_report_coordinates} = 'true';

    my $annotation_info->{data_format} = 'vcf';
    $annotation_info->{file} = { class => 'File', path => $vcf_path, secondaryFiles => [{class => 'File', path => $vcf_path . '.tbi' }]};
    $annotation_info->{name} = 'clinvar';
    $annotation_info->{gnomad_filter} = 'false';
    $annotation_info->{check_existing} = 'false';
    $annotation_info->{vcf_fields} = ['CLINSIGN','PHENOTYPE','SCORE','RCVACC','TESTEDINGTR','PHENOTYPELIST','NUMSUBMIT','GUIDELINES'];
    $custom_annotation->{annotation} = $annotation_info;

    push @vep_custom_annotations, $custom_annotation;
}
$inputs->{vep_custom_annotations} = \@vep_custom_annotations;

my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
YAML::XS::DumpFile($yaml, $inputs);
