#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: "STAR: align reads to transcriptome"
baseCommand: ["/usr/local/bin/STAR"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 12
    - class: DockerRequirement
      dockerPull: "trinityctat/starfusion:1.10.1"
arguments: [
    "--runThreadN", $(runtime.cores)

]
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
    out_reads_unmapped:
        type: string
        default: "None"
        inputBinding:
            position: 6
            prefix: '--outReadsUnmapped'
    chim_segment_min:
        type: int
        default: "12"
        inputBinding:
            position: 7
            prefix: '--chimSegmentMin'
    chim_junction_overhang_min:
        type: int
        default: "12"
        inputBinding:
            position: 8
            prefix: '--chimJunctionOverhangMin'
    align_sjdb_overhang_min:
        type: int
        default: 10
        inputBinding:
            position: 9
            prefix: '--alignSJDBoverhangMin'
    align_mates_gapmax:
        type: int?
        default: 100000
        inputBinding:
            position: 10
            prefix: '--alignMatesGapMax'
    align_intron_max:
        type: int
        default: 100000
        inputBinding:
            position: 11
            prefix: '--alignIntronMax'
    chim_segment_read_gapmax:
        type: int
        default: 3
        inputBinding:
            position: 12
            prefix: '--chimSegmentReadGapMax'
    align_sjstitch_mismatch_nmax:
        type: int[]
        default: [5, -1, 5, 5]
        inputBinding:
            position: 13 
            prefix: '--alignSJstitchMismatchNmax'
            itemSeparator: ' '
            shellQuote: False
    outsam_strand_field:
        type: string
        default: "intronMotif"
        inputBinding:
            position: 14
            prefix: '--outSAMstrandField'
    outsam_unmapped:
        type: string
        default: Within
        inputBinding:
            position: 15
            prefix: '--outSAMunmapped'
    outsam_attrrg_line:
        type:
            type: array
            items: string
        inputBinding:
            position: 16
            itemSeparator: ' , ' 
            shellQuote: False
            prefix: '--outSAMattrRGline'
        doc: '
            string(s): SAM/BAM read group line. The first word contains the read group
            identifier and must start with ID:, e.g. --outSAMattrRGline ID:xxx CN:yy
            DS:z z z.
            xxx will be added as RG tag to each output alignment. Any spaces in the tag
            values have to be double quoted.
            Comma separated RG lines correspons to different (comma separated) input
            files in --readFilesIn
            '
    chim_multimap_nmax:
        type: int
        default: 10
        inputBinding:
            position: 17
            prefix: '--chimMultimapNmax'
    chim_nonchim_scoredrop_min:
        type: int
        default: 10
        inputBinding:
            position: 18
            prefix: '--chimNonchimScoreDropMin'
    peoverlap_nbases_min:
        type: int
        default: 12
        inputBinding:
            position: 19
            prefix: '--peOverlapNbasesMin'
    peoverlap_mmp:
        type: float
        default: 0.1
        inputBinding:
            position: 20
            prefix: '--peOverlapMMp'
    chimout_junction_format:
        type: int
        default: 1
        inputBinding:
            position: 21
            prefix: '--chimOutJunctionFormat'
    star_genome_dir:
        type: Directory
        inputBinding:
            position: 22
            prefix: '--genomeDir'
        doc: '
            specifies path to the directory where the genome indices are stored
            '
    twopass_mode:
        type: string
        default: "Basic"
        inputBinding:
            position: 23
            prefix: '--twopassMode'
    reference_annotation:
        type: File?
        inputBinding:
            position: 24
            prefix: '--sjdbGTFfile'
        doc: 'Annotated transcripts in GTF format; used as a source of splice junctions'
    outfile_name_prefix:
        type: string
        default: "STAR_"
        inputBinding:
            position: 25
            prefix: '--outFileNamePrefix'
    read_files_command:
        default: "cat"
        type: string?
        inputBinding:
            position: 26
            prefix: '--readFilesCommand'
    outsam_attributes:
        type: string[]
        default: [NH, HI, AS, NM, MD]
        inputBinding:
            position: 27
            prefix: '--outSAMattributes'
            itemSeparator: ' '
            shellQuote: False

outputs:
    aligned_bam:
        type: File
        outputBinding:
          glob: "$(inputs.outfile_name_prefix)Aligned.out.bam"
    log_final:
        type: File
        outputBinding:
          glob: "$(inputs.outfile_name_prefix)Log.final.out"
    log:
        type: File
        outputBinding:
          glob: "$(inputs.outfile_name_prefix)Log.out"
    log_progress:
        type: File
        outputBinding:
          glob: "$(inputs.outfile_name_prefix)Log.progress.out"
    splice_junction_out:
        type: File
        outputBinding:
            glob: "$(inputs.outfile_name_prefix)SJ.out.tab"
    chim_junc:
        type: File
        outputBinding:
            glob: "$(inputs.outfile_name_prefix)Chimeric.out.junction"
