#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: "STAR: align reads to transcriptome for Xenosplit"
baseCommand: ["/usr/local/bin/STAR"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 10
    - class: DockerRequirement
      dockerPull: "mgibio/star:2.7.0f"
arguments: ["--runThreadN", $(runtime.cores)]

inputs:
    out_samtype:
        type: string[]
        default: ["BAM", "Unsorted"]
        inputBinding:
            position: 3
            prefix: '--outSAMtype'
    run_mode:
        type: string
        default: "alignReads"
        inputBinding:
            position: 2
            prefix: "--runMode"
    fastq:
        type: File[]
        inputBinding:
            position: 4
            prefix: '--readFilesIn'
            itemSeparator: ","
    fastq2:
          type: File[]
          inputBinding:
              position: 5
              itemSeparator: ","
    outsam_unmapped:
        type: string
        default: Within
        inputBinding:
            position: 6
            prefix: '--outSAMunmapped'
    star_genome_dir:
        type: Directory
        inputBinding:
            position: 7
            prefix: '--genomeDir'
        doc: '
            specifies path to the directory where the genome indices are stored
            '
    outfile_name_prefix:
        type: string
        inputBinding:
            position: 8
            prefix: '--outFileNamePrefix'

outputs:
    aligned_bam:
        type: File
        outputBinding:
          glob: "$(inputs.outfile_name_prefix)Aligned.out.bam"
