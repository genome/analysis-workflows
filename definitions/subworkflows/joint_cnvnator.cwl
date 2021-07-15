#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "run cnvnator for multiple samples"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    sample_names:
        type: string[]
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    bin_size:
        type: int?
outputs:
    vcfs:
        type: File[]
        outputSource: index_cnvnator/indexed_vcf
        secondaryFiles: [.tbi]
    root_files:
        type: File[]
        outputSource: cnvnator/root_file
    cn_files:
        type: File[]
        outputSource: cnvnator/cn_file
steps:
    cnvnator:
        scatter: [bam, sample_name]
        scatterMethod: dotproduct
        run: ../tools/cnvnator.cwl
        in:
            bam: bams
            reference: reference
            sample_name: sample_names
            bin_size: bin_size
        out:
            [vcf, root_file, cn_file]
    bgzip_index:
        scatter: [vcf]
        run: bgzip_and_index.cwl
        in:
            vcf: cnvnator/vcf
        out:
            [indexed_vcf]
    sample_rename:
        scatter: [input_vcf, new_sample_name]
        scatterMethod: dotproduct
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: bgzip_index/indexed_vcf
            new_sample_name: sample_names
            sample_to_replace:
                valueFrom: '${
                    var old_name = inputs.new_sample_name.split(".")[0];
                    return old_name;
                }'
            output_name:
                valueFrom: '${
                    var sample = inputs.new_sample_name;
                    var name = sample + ".cnvnator.vcf.gz";
                    return name;
                }'
        out:
            [renamed_vcf]
    index_cnvnator:
        scatter: [vcf]
        run: ../tools/index_vcf.cwl
        in:
            vcf: sample_rename/renamed_vcf
        out:
            [indexed_vcf]
