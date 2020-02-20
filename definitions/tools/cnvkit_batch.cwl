#! /usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["/usr/bin/python", "/usr/local/bin/cnvkit.py", "batch"]

requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "etal/cnvkit:0.9.5"
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 4000
      tmpdirMin: 10000
    - class: InlineJavascriptRequirement
arguments: [{ valueFrom: "$((inputs.reference.hasOwnProperty('cnn_file'))? null : '--normal')" }]
inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: -1
    bait_intervals:
        type: File?
        inputBinding:
            position: 3
            prefix: "--targets"
    reference:
      type:
        - type: record
          name: cnn_file
          fields:
            cnn_file:
              type: File
              inputBinding:
                position: 2
                prefix: "--reference"
              doc: "Previously generated reference.cnn file"
        - type: record
          name: fasta_file
          fields:
            fasta_file:
              type:
                - string
                - File
              inputBinding:
                position: 2
                prefix: "--fasta"
            normal_bam:
              type: File?
              inputBinding:
                position: 1
              doc: "Normal samples (.bam) used to construct the pooled, paired, or flat reference. If this option is used but no filenames are given, a 'flat' reference will be built. Otherwise, all filenames following this option will be used."
    access:
        type: File?
        inputBinding:
            position: 4
            prefix: "--access"
        doc: "Regions of accessible sequence on chromosomes (.bed), as output by the 'access' command"
    method:
        type:
            - "null"
            - type: enum
              symbols: ["hybrid", "amplicon", "wgs"]
        default: "hybrid"
        inputBinding:
            position: 5
            prefix: "--method"
        doc: "Sequencing protocol used for input data"
    diagram:
        type: boolean?
        inputBinding:
            position: 6
            prefix: "--diagram"
        doc: "Create an ideogram of copy ratios on chromosomes as a PDF"
    scatter_plot:
        type: boolean?
        inputBinding:
            position: 7
            prefix: "--scatter"
        doc: "Create a whole-genome copy ratio profile as a PDF scatter plot"
    drop_low_coverage:
        type: boolean?
        inputBinding:
            position: 8
            prefix: "--drop-low-coverage"
        doc: "Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples"
    male_reference:
        type: boolean?
        inputBinding:
            position: 9
            prefix: "--male-reference"
        doc: "Use or assume a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX)"
    target_average_size:
        type: int?
        inputBinding:
            position: 10
            prefix: "--target-avg-size"
        doc: "Average size of split target bins (results are approximate)"
outputs:
    intervals_antitarget:
        type: File?
        outputBinding:
            glob: |
                    ${  
                        var glob_base = ".antitarget.bed";
                        if (inputs.bait_intervals) {
                            glob_base = inputs.bait_intervals.nameroot + glob_base;
                        }   
                        return glob_base;
                    }  
    intervals_target:
        type: File?
        outputBinding:
            glob: |
                    ${
                        var glob_base = ".target.bed";
                        if (inputs.bait_intervals) {
                            glob_base = inputs.bait_intervals.nameroot + glob_base;
                        }
                        return glob_base;
                    }
    normal_antitarget_coverage:
        type: File?
        outputBinding:
            glob: |
                    ${
                        var glob_base = ".antitargetcoverage.cnn";
                        if (inputs.normal_bam) {
                            glob_base = inputs.normal_bam.nameroot + glob_base;
                        }
                        return glob_base;
                    }
    normal_target_coverage:
        type: File?
        outputBinding:
            glob: |
                    ${
                        var glob_base = ".targetcoverage.cnn";
                        if (inputs.normal_bam) {
                            glob_base = inputs.normal_bam.nameroot + glob_base;
                        }
                        return glob_base;
                    }
    reference_coverage:
        type: File?
        outputBinding:
            glob: reference.cnn
    cn_diagram:
        type: File?
        outputBinding:
            glob: $(inputs.tumor_bam.nameroot)-diagram.pdf
    cn_scatter_plot:
        type: File?
        outputBinding:
            glob: $(inputs.tumor_bam.nameroot)-scatter.pdf
    tumor_antitarget_coverage:
        type: File
        outputBinding:
            glob: $(inputs.tumor_bam.nameroot).antitargetcoverage.cnn
    tumor_target_coverage:
        type: File
        outputBinding:
            glob: $(inputs.tumor_bam.nameroot).targetcoverage.cnn
    tumor_bin_level_ratios:
        type: File
        outputBinding:
            glob: $(inputs.tumor_bam.nameroot).cnr
    tumor_segmented_ratios:
        type: File
        outputBinding:
            glob: $(inputs.tumor_bam.nameroot).cns
