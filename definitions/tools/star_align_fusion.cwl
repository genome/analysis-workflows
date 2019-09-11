#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: "STAR: align reads to transcriptome"
baseCommand: ["/bin/bash star.sh"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 10
    - class: DockerRequirement
      dockerPull: "mgibio/star:2.7.0f"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'star.sh'
        entry: |
            
            set -eou pipefail

            passthrough_args=${@:1:25}
            paired="$26"
            fastqs="$27"
            
            if [[ "$paired" == "false" ]];then
                #run star with all the fastqs in single-end mode
                /usr/local/bin/STAR "$passthrough_args" --readFilesIn $(join , $fastqs)
            else 
                #split interleaved fastqs into fq1/2 arguments
                i=0
                fastq_array=(${fastqs})
                fq1=()
                fq2=()
                while [[ "$i" -lt ${#fastq_array[*]} ]];do
                    fq1+=( "${fastq_array[$i]}" )
                    fq2+=( "${fastq_array[$(($i+1))]}" )
                    let "i=$i+2"
                done
                /usr/local/bin/STAR "$passthrough_args" --readFilesIn $(join , "${fq1[@]}") $(join , "${fq2[@]}")

arguments: [
    $(runtime.cores)
]
inputs:
    run_mode:
        type: string
        default: "alignReads"
        inputBinding:
            position: 2
            prefix: "--runMode"
    out_samtype:
        type: string[]
        default: ["BAM", "Unsorted"]
        inputBinding:
            position: 3
            prefix: '--outSAMtype'
    out_reads_unmapped:
        type: string
        default: "None"
        inputBinding:
            position: 4
            prefix: '--outReadsUnmapped'
    chim_segment_min:
        type: int
        default: "12"
        inputBinding:
            position: 5
            prefix: '--chimSegmentMin'
    chim_junction_overhang_min:
        type: int
        default: "12"
        inputBinding:
            position: 6
            prefix: '--chimJunctionOverhangMin'
    align_sjdb_overhang_min:
        type: int
        default: 10
        inputBinding:
            position: 7
            prefix: '--alignSJDBoverhangMin'
    align_mates_gapmax:
        type: int?
        default: 100000
        inputBinding:
            position: 8
            prefix: '--alignMatesGapMax'
    align_intron_max:
        type: int
        default: 100000
        inputBinding:
            position: 9
            prefix: '--alignIntronMax'
    chim_segment_read_gapmax:
        type: int
        default: 3
        inputBinding:
            position: 10
            prefix: '--chimSegmentReadGapMax'
    align_sjstitch_mismatch_nmax:
        type: int[]
        default: [5, -1, 5, 5]
        inputBinding:
            position: 11
            prefix: '--alignSJstitchMismatchNmax'
            itemSeparator: ' '
            shellQuote: False
    outsam_strand_field:
        type: string
        default: "intronMotif"
        inputBinding:
            position: 12
            prefix: '--outSAMstrandField'
    outsam_unmapped:
        type: string
        default: Within
        inputBinding:
            position: 13
            prefix: '--outSAMunmapped'
    outsam_attrrg_line:
        type:
            type: array
            items: string
        inputBinding:
            position: 14
            itemSeparator: ' , ' 
            shellQuote: False
            prefix: '--outSAMattrRGline'
        doc: '
            string(s): SAM/BAM read group line. The first word contains the read group
            identifier and must start with ID:, e.g. –outSAMattrRGline ID:xxx CN:yy
            DS:z z z.
            xxx will be added as RG tag to each output alignment. Any spaces in the tag
            values have to be double quoted.
            Comma separated RG lines correspons to different (comma separated) input
            files in –readFilesIn
            '
    chim_multimap_nmax:
        type: int
        default: 10
        inputBinding:
            position: 15
            prefix: '--chimMultimapNmax'
    chim_nonchim_scoredrop_min:
        type: int
        default: 10
        inputBinding:
            position: 16
            prefix: '--chimNonchimScoreDropMin'
    peoverlap_nbases_min:
        type: int
        default: 12
        inputBinding:
            position: 17
            prefix: '--peOverlapNbasesMin'
    peoverlap_mmp:
        type: float
        default: 0.1
        inputBinding:
            position: 18
            prefix: '--peOverlapMMp'
    chimout_junction_format:
        type: int
        default: 1
        inputBinding:
            position: 19
            prefix: '--chimOutJunctionFormat'
    star_genome_dir:
        type: Directory
        inputBinding:
            position: 20
            prefix: '--genomeDir'
        doc: '
            specifies path to the directory where the genome indices are stored
            '
    twopass_mode:
        type: string
        default: "Basic"
        inputBinding:
            position: 21
            prefix: '--twopassMode'
    gtf_file:
        type: File?
        inputBinding:
            position: 22
            prefix: '--sjdbGTFfile'
    outfile_name_prefix:
        type: string
        default: "STAR_"
        inputBinding:
            position: 23
            prefix: '--outFileNamePrefix'
    read_files_command:
        default: "cat"
        type: string?
        inputBinding:
            position: 24
            prefix: '--readFilesCommand'
    outsam_attributes:
        type: string[]
        default: [NH, HI, AS, NM, MD]
        inputBinding:
            position: 25
            prefix: '--outSAMattributes'
            itemSeparator: ' '
            shellQuote: False
    paired_end:
        type: boolean
        default: true
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
        inputBinding:
            position: 26
    fastqs:
        type: File[]
        inputBinding:
            position: 27
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
