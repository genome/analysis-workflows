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
inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: 1
    normal_bam:
        type: File?
        inputBinding:
            position: 2
            prefix: "--normal"
    bait_intervals:
        type: File?
        inputBinding:
            position: 3
            prefix: "--targets"
    reference:
        type: string?
        inputBinding:
            position: 4
            prefix: "--fasta"
    access:
        type: File?
        inputBinding:
            position: 5
            prefix: "--access"
    method:
        type: string?
        default: "hybrid"
        inputBinding:
            position: 6
            prefix: "--method"
    diagram:
        type: boolean?
        inputBinding:
            position: 7
            prefix: "--diagram"
    scatter_plot:
        type: boolean?
        inputBinding:
            position: 8
            prefix: "--scatter"
    drop_low_coverage:
        type: boolean?
        inputBinding:
            position: 10
            prefix: "--drop-low-coverage"
    male_reference:
        type: boolean?
        inputBinding:
            prefix: "-y"
    reference_cnn:
        type: File?
        inputBinding:
            prefix: "-r"
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
