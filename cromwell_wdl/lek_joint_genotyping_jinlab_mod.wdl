## Copyright Broad Institute, 2018
## 
## This WDL implements the joint discovery and VQSR filtering portion of the GATK 
## Best Practices (June 2016) for germline SNP and Indel discovery in human 
## whole-genome sequencing (WGS) and exome sequencing data.
##
## Requirements/expectations :
## - One or more GVCFs produced by HaplotypeCaller in GVCF mode 
## - Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
##
## Outputs :
## - A VCF file and its index, filtered using variant quality score recalibration 
##   (VQSR) with genotypes for all samples present in the input VCF. All sites that 
##   are present in the input VCF are retained; filtered sites are annotated as such 
##   in the FILTER field.
##
## Note about VQSR wiring :
## The SNP and INDEL models are built in parallel, but then the corresponding 
## recalibrations are applied in series. Because the INDEL model is generally ready 
## first (because there are fewer indels than SNPs) we set INDEL recalibration to 
## be applied first to the input VCF, while the SNP model is still being built. By 
## the time the SNP model is available, the indel-recalibrated file is available to 
## serve as input to apply the SNP recalibration. If we did it the other way around, 
## we would have to wait until the SNP recal file was available despite the INDEL 
## recal file being there already, then apply SNP recalibration, then apply INDEL 
## recalibration. This would lead to a longer wall clock time for complete workflow 
## execution. Wiring the INDEL recalibration to be applied first solves the problem.
##
## Cromwell version support 
## - Successfully tested on v31
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

## Adapted to Yale Ruddle HPC by Sander Pajusalu (sander.pajusalu@yale.edu)
## Modified for Washington-University compute0 and compute1 by Samuel Peters 6/19/20
## Modifed to use GATK 4.1.7.0 and Docker

## DOCKER
## Referenced digest ID for exact container instead of tag
## broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2

workflow JointGenotyping {
  String unpadded_intervals_file

  String callset_name
  
  String ref_fasta
  String ref_fasta_index
  String ref_dict

  String dbsnp_vcf
  String dbsnp_vcf_index

  String sample_sheet

  Array[String] snp_recalibration_tranche_values
  Array[String] snp_recalibration_annotation_values
  Array[String] indel_recalibration_tranche_values
  Array[String] indel_recalibration_annotation_values

  String eval_interval_list
  String hapmap_resource_vcf
  String hapmap_resource_vcf_index
  String omni_resource_vcf
  String omni_resource_vcf_index
  String one_thousand_genomes_resource_vcf
  String one_thousand_genomes_resource_vcf_index
  String mills_resource_vcf
  String mills_resource_vcf_index
  String axiomPoly_resource_vcf
  String axiomPoly_resource_vcf_index
  String dbsnp_resource_vcf = dbsnp_vcf
  String dbsnp_resource_vcf_index = dbsnp_vcf_index

  String DynamicallyCombineIntervals_py
  String statsmerge_v2_washu

  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69
  Float snp_filter_level
  Float indel_filter_level
  Int SNP_VQSR_downsampleFactor

  Int num_of_original_intervals = length(read_lines(unpadded_intervals_file))
  

  # Make a 2.5:1 interval number to samples in callset ratio interval list
  Int possible_merge_count = floor(num_of_original_intervals / num_gvcfs / 2.5)
  Int merge_count = if possible_merge_count > 1 then possible_merge_count else 1

  call samples {
    input:
      samples = sample_sheet
  }
  
  call exomeMetrics {
    input:
      statsmerge_v2_washu = statsmerge_v2_washu,
      metrics_paths = samples.results_paths
  }

  Int num_gvcfs = length(read_lines(samples.input_gvcfs))

  call DynamicallyCombineIntervals {
    input:
      intervals = unpadded_intervals_file,
      merge_count = merge_count,
      dyn_combine_int = DynamicallyCombineIntervals_py

  }

  Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)

  scatter (idx in range(length(unpadded_intervals))) {
    # the batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call ImportGVCFs {
      input:
        sample_names = read_lines(samples.sample_names),
        interval = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        input_gvcfs = read_lines(samples.input_gvcfs),
        input_gvcfs_indices = read_lines(samples.input_gvcfs_indices),
        batch_size = 50
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index
    }

    call HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        vcf_index = GenotypeGVCFs.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz"
    }
  }

  call GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs_fofn = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
      input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index,
      output_vcf_name = callset_name + ".sites_only.vcf.gz"
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
  }

  if (num_gvcfs > 10000) {
  call SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        downsampleFactor = SNP_VQSR_downsampleFactor,
        model_report_filename = callset_name + ".snps.model.report",
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
    }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
    call SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
      input:
        sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
        sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
        recalibration_filename = callset_name + ".snps." + idx + ".recal",
        tranches_filename = callset_name + ".snps." + idx + ".tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        model_report = SNPsVariantRecalibratorCreateModel.model_report,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
      }
    }
    call GatherTranches as SNPGatherTranches {
        input:
          input_fofn = SNPsVariantRecalibratorScattered.tranches,
          output_filename = callset_name + ".snps.gathered.tranches"
    }
  }

  if (num_gvcfs <= 10000){
    call SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
          sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
          sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
          recalibration_filename = callset_name + ".snps.recal",
          tranches_filename = callset_name + ".snps.tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
    }
  }

  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  Boolean is_small_callset = num_gvcfs <= 1000

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    call ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
        snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
        snps_tranches = select_first([SNPGatherTranches.tranches, SNPsVariantRecalibratorClassic.tranches]),
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level
    }

    # for large callsets we need to collect metrics from the shards and gather them later
    if (!is_small_callset) {
      call CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict
      }
    }
  }

  # for small callsets we can gather the VCF shards and then collect metrics on it
  if (is_small_callset) {
    call GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs_fofn = ApplyRecalibration.recalibrated_vcf,
        input_vcf_indexes_fofn = ApplyRecalibration.recalibrated_vcf_index,
        output_vcf_name = callset_name + ".vcf.gz"
    }

    call CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict
    }
  }

  # for large callsets we still need to gather the sharded metrics
  if (!is_small_callset) {
    call GatherMetrics {
      input:
        input_details_fofn = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries_fofn = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = callset_name
    }
  }

  output {
    # outputs exome metrics summary file for bamMetrics results from each sample
    exomeMetrics.exome_metrics
    # outputs from the small callset path through the wdl
    FinalGatherVcf.output_vcf
    FinalGatherVcf.output_vcf_index
    CollectMetricsOnFullVcf.detail_metrics_file
    CollectMetricsOnFullVcf.summary_metrics_file

    # outputs from the large callset path through the wdl
    # (note that we do not list ApplyRecalibration here because it is run in both paths)
    GatherMetrics.detail_metrics_file
    GatherMetrics.summary_metrics_file

    # output the interval list generated/used by this run workflow
    DynamicallyCombineIntervals.output_intervals
  }
}

task samples {
  File samples

  command {
    cut ${samples} -f1 > sample_names.txt
    cut ${samples} -f2 > gvcfs.txt
    cut ${samples} -f3 > gvcf_indices.txt
    cut ${samples} -f4 > results_paths.txt

  }
  runtime {
    cpus: 2
    requested_memory: 4000  
  }
  output {
    File sample_names = "sample_names.txt"
    File input_gvcfs = "gvcfs.txt"
    File input_gvcfs_indices = "gvcf_indices.txt"
    File results_paths = "results_paths.txt"
  }

}

task exomeMetrics{
  File metrics_paths
  String statsmerge_v2_washu
  command {
      python ${statsmerge_v2_washu} \
      ${metrics_paths} > exome_metrics.txt
  }
  runtime {
    cpus: 2
    requested_memory: 4000  
  }
  output {
    File exome_metrics = stdout()
  }
}

task GetNumberOfSamples {
  File sample_name_map
  command <<<
    wc -l ${sample_name_map} | awk '{print $1}'
  >>>

  runtime {
    cpus: 4
    requested_memory: 8000
  }
  output {
    Int sample_count = read_int(stdout())
  }
}

task ImportGVCFs {
  Array[String] sample_names
  Array[String] input_gvcfs
  Array[String] input_gvcfs_indices
  String interval

  String workspace_dir_name

  Int batch_size

  command <<<
    set -e
    set -o pipefail
    
    python << CODE
    gvcfs = ['${sep="','" input_gvcfs}']
    sample_names = ['${sep="','" sample_names}']

    if len(gvcfs)!= len(sample_names):
      exit(1)

    with open("inputs.list", "w") as fi:
      for i in range(len(gvcfs)):
        fi.write(sample_names[i] + "\t" + gvcfs[i] + "\n") 

    CODE

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    /gatk/gatk --java-options "-Xmx4g -Xms4g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${interval} \
    --sample-name-map inputs.list \
    --reader-threads 5 \
    -ip 500

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

  >>>
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  String ref_fasta
  String ref_fasta_index
  String ref_dict

  String dbsnp_vcf
  String dbsnp_vcf_index


  command <<<
    set -e

    tar -xf ${workspace_tar}
    WORKSPACE=$( basename ${workspace_tar} .tar)

    /gatk/gatk --java-options "-Xmx14g -Xms5g" \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbsnp_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     -V gendb://$WORKSPACE \
     -L ${interval}
  >>>
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 16000  
  }
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename


  command {
    set -e

    /gatk/gatk --java-options "-Xmx15g -Xms3g" \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    /gatk/gatk --java-options "-Xmx15g -Xms3g" \
      MakeSitesOnlyVcf \
      --INPUT ${variant_filtered_vcf_filename} \
      --OUTPUT ${sites_only_vcf_filename}

  }
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 16000  
  }
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  }
}

task IndelsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  String mills_resource_vcf
  String axiomPoly_resource_vcf
  String dbsnp_resource_vcf
  String mills_resource_vcf_index
  String axiomPoly_resource_vcf_index
  String dbsnp_resource_vcf_index


  command {
    /gatk/gatk --java-options "-Xmx24g -Xms24g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode INDEL \
      --max-gaussians 4 \
      --resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
      --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf}
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 8
    requested_memory: 32000
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task SNPsVariantRecalibratorCreateModel {
  String recalibration_filename
  String tranches_filename
  Int downsampleFactor
  String model_report_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  String hapmap_resource_vcf
  String omni_resource_vcf
  String one_thousand_genomes_resource_vcf
  String dbsnp_resource_vcf
  String hapmap_resource_vcf_index
  String omni_resource_vcf_index
  String one_thousand_genomes_resource_vcf_index
  String dbsnp_resource_vcf_index


  command {
    /gatk/gatk --java-options "-Xmx100g -Xms100g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      --sample-every-Nth-variant ${downsampleFactor} \
      --output-model ${model_report_filename} \
      --max-gaussians 6 \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
      --resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 8
    requested_memory: 104000
  }
  output {
    File model_report = "${model_report_filename}"
  }
}

task SNPsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename
  File? model_report

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  String hapmap_resource_vcf
  String omni_resource_vcf
  String one_thousand_genomes_resource_vcf
  String dbsnp_resource_vcf
  String hapmap_resource_vcf_index
  String omni_resource_vcf_index
  String one_thousand_genomes_resource_vcf_index
  String dbsnp_resource_vcf_index


  command {
    /gatk/gatk --java-options "-Xmx3g -Xms3g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      ${"--input-model " + model_report + " --output-tranches-for-scatter "} \
      --max-gaussians 6 \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
      --resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}

task GatherTranches {
  Array[File] input_fofn
  String output_filename


  command <<<
    set -e
    set -o pipefail
    
      /gatk/gatk --java-options "-Xmx6g -Xms6g" \
      GatherTranches \
      --input ${sep=" --input " input_fofn}  \
      --output ${output_filename}
  >>>
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
  output {
    File tranches = "${output_filename}"
  }
}

task ApplyRecalibration {
  String recalibrated_vcf_filename
  String input_vcf
  String input_vcf_index
  File indels_recalibration
  File indels_recalibration_index
  File indels_tranches
  File snps_recalibration
  File snps_recalibration_index
  File snps_tranches

  Float indel_filter_level
  Float snp_filter_level

  command {
    set -e

    /gatk/gatk --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    /gatk/gatk --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O ${recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snps_recalibration} \
      --tranches-file ${snps_tranches} \
      --truth-sensitivity-filter-level ${snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
  output {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}

task GatherVcfs {
  Array[File] input_vcfs_fofn
  Array[File] input_vcf_indexes_fofn
  String output_vcf_name
  
  command <<<
    set -e
    set -o pipefail

    # ignoreSafetyChecks make a big performance difference so we include it in our invocation
    /gatk/gatk --java-options "-Xmx6g -Xms6g" \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    --input ${sep=" --input " input_vcfs_fofn} \
    --output ${output_vcf_name}

    /gatk/gatk --java-options "-Xmx6g -Xms6g" \
    IndexFeatureFile \
    --input ${output_vcf_name}
  >>>
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {
  String input_vcf
  String input_vcf_index

  String metrics_filename_prefix
  String dbsnp_vcf
  String dbsnp_vcf_index
  String interval_list
  String ref_dict


  command {
    /gatk/gatk --java-options "-Xmx6g -Xms6g" \
      CollectVariantCallingMetrics \
      --INPUT ${input_vcf} \
      --DBSNP ${dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ${ref_dict} \
      --OUTPUT ${metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${interval_list}
  }
  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
}

task GatherMetrics {
  Array[File] input_details_fofn
  Array[File] input_summaries_fofn

  String output_prefix

  command <<<
    set -e
    set -o pipefail

    
    /gatk/gatk --java-options "-Xmx2g -Xms2g" \
    AccumulateVariantCallingMetrics \
    --INPUT ${sep=" --INPUT " input_details_fofn} \
    --OUTPUT ${output_prefix}
  >>>
  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 8000
  }
  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}

task DynamicallyCombineIntervals {
  File intervals
  Int merge_count
  File dyn_combine_int

  #command {
    #python << CODE
    #def parse_interval(interval):
        #colon_split = interval.split(":")
        #chromosome = colon_split[0]
        #dash_split = colon_split[1].split("-")
        #start = int(dash_split[0])
        #end = int(dash_split[1])
        #return chromosome, start, end

    #def add_interval(chr, start, end):
        #lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        #return chr, start, end

    #count = 0
    #chain_count = ${merge_count}
    #l_chr, l_start, l_end = "", 0, 0
    #lines_to_write = []
    #with open("${intervals}") as f:
        #with open("out.intervals", "w") as f1:
            #for line in f.readlines():
                # initialization
                #if count == 0:
                    #w_chr, w_start, w_end = parse_interval(line)
                    #count = 1
                    #continue
                # reached number to combine, so spit out and start over
                #if count == chain_count:
                    #l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    #w_chr, w_start, w_end = parse_interval(line)
                    #count = 1
                    #continue

                #c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                #if c_chr == w_chr and c_start == w_end + 1:
                    #w_end = c_end
                    #count += 1
                    #continue
                # not adjacent, end here and start a new chain
                #else:
                    #l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    #w_chr, w_start, w_end = parse_interval(line)
                    #count = 1
            #if l_char != w_chr or l_start != w_start or l_end != w_end:
                #add_interval(w_chr, w_start, w_end)
            #f1.writelines("\n".join(lines_to_write))
    #CODE
  #}
  command {
      python ${dyn_combine_int} \
      ${intervals} \
      ${merge_count}
  }

  runtime {
    docker: "broadinstitute/gatk:4.1.7.0@sha256:192fedf9b9d65809da4a2954030269e3f311d296e6b5e4c6c7deec12c7fe84b2"
    cpus: 4
    requested_memory: 15000
  }

  output {
    File output_intervals = "out.intervals"
  }
}