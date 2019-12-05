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

            while getopts "c:r:i:g:p:s:f:" opt; do
                case "$opt" in
                    c)
                        cores="$OPTARG"
                        ;;
                case "$opt" in
                    f)
                        fastqs="$OPTARG"
                        ;;
                case "$opt" in
                    s)
                        chimSegmentMin="$OPTARG"
                        ;;
                case "$opt" in
                    j)
                        chimJunctionOverhangMin="$OPTARG"
                        ;;
                case "$opt" in
                    b)
                        alignSJDBoverhangMin="$OPTARG"
                        ;;
                case "$opt" in
                    q)
                        alignMatesGapMax="$OPTARG"
                        ;;
                case "$opt" in
                    i)
                        alignIntronMax="$OPTARG"
                        ;;
                case "$opt" in
                    r)
                        chimSegmentReadGapMax="$OPTARG"
                        ;;
                case "$opt" in
                    m)
                        alignSJstitchMismatchNmax="$OPTARG"
                        ;;
                case "$opt" in
                    t)
                        outSAMstrandField="$OPTARG"
                        ;;
                case "$opt" in
                    u)
                        outSAMunmapped="$OPTARG"
                        ;;
                case "$opt" in
                    a)
                        outSAMattrRGline="$OPTARG"
                        ;;
                case "$opt" in
                    h)
                        chimMultimapNmax="$OPTARG"
                        ;;
                case "$opt" in
                    d)
                        chimNonchimScoreDropMin="$OPTARG"
                        ;;
                case "$opt" in
                    o)
                        peOverlapNbasesMin="$OPTARG"
                        ;;
                case "$opt" in
                    p)
                        peOverlapMMp="$OPTARG"
                        ;;
                case "$opt" in
                    v)
                        chimOutJunctionFormat="$OPTARG"
                        ;;
                case "$opt" in
                    g)
                        genomeDir="$OPTARG"
                        ;;
                case "$opt" in
                    w)
                        twopassMode="$OPTARG"
                        ;;
                case "$opt" in
                    k)
                        sjdbGTFfile="$OPTARG"
                        ;;
                case "$opt" in
                    e)
                        paired="$OPTARG"
                        ;;

                esac
            done

            fqfinal=""
            if [[ "$paired" == "false" ]];then
                #run star with all the fastqs in single-end mode
                fqfinal=`join , $fastqs`
            else
                #split interleaved fastqs into fq1/2 arguments
                i=0
                fastq_array=(${fastqs})
                fq1=()
                fq2=()
                while [[ "$i" -lt ${#fastq_array[*]} ]];do
                    fq1+=( "${fastq_array[$i]}" )
                    let "i=$i+1"
                    fq2+=( "${fastq_array[$i]}" )
                    let "i=$i+1"
                done
                fqfinal=`join , "${fq1[@]}"` `join , "${fq2[@]}"`
            fi

            /usr/local/bin/STAR --runMode alignReads --outSAMtype BAM Unsorted --outReadsUnmapped None --outFileNamePrefix STAR_ --readFilesCommand cat --outSAMattributes NH HI AS NM MD --runThreadN "$cores" --readFilesIn "$fqfinal"--chimSegmentMin "$chimSegmentMin" --chimJunctionOverhangMin "$chimJunctionOverhangMin" --alignSJDBoverhangMin "$alignSJDBoverhangMin" --alignMatesGapMax "$alignMatesGapMax" --alignIntronMax "$alignIntronMax" --chimSegmentReadGapMax "$chimSegmentReadGapMax" --alignSJstitchMismatchNmax "$alignSJstitchMismatchNmax" --outSAMstrandField "$outSAMstrandField" --outSAMunmapped "$outSAMunmapped" --outSAMattrRGline "$outSAMattrRGline" --chimMultimapNmax "$chimMultimapNmax" --chimNonchimScoreDropMin "$chimNonchimScoreDropMin" --peOverlapNbasesMin "$peOverlapNbasesMin" --peOverlapMMp "$peOverlapMMp" --chimOutJunctionFormat "$chimOutJunctionFormat" --genomeDir "$genomeDir" --twopassMode "$twopassMode" --sjdbGTFfile "$sjdbGTFfile"

arguments: [
    {valueFrom: "$(runtime.cores)", position: 1, prefix: "-c"}
]
inputs:
    chim_segment_min:
        type: int
        default: "12"
        inputBinding:
            position: 5
            prefix: '-s'
    chim_junction_overhang_min:
        type: int
        default: "12"
        inputBinding:
            position: 6
            prefix: '-j'
    align_sjdb_overhang_min:
        type: int
        default: 10
        inputBinding:
            position: 7
            prefix: '-b'
    align_mates_gapmax:
        type: int?
        default: 100000
        inputBinding:
            position: 8
            prefix: '-q'
    align_intron_max:
        type: int
        default: 100000
        inputBinding:
            position: 9
            prefix: '-i'
    chim_segment_read_gapmax:
        type: int
        default: 3
        inputBinding:
            position: 10
            prefix: '-r'
    align_sjstitch_mismatch_nmax:
        type: int[]
        default: [5, -1, 5, 5]
        inputBinding:
            position: 11
            prefix: '-m'
            itemSeparator: ' '
            shellQuote: False
    outsam_unmapped:
        type: string
        default: Within
        inputBinding:
            position: 13
            prefix: '-u'
    outsam_attrrg_line:
        type:
            type: array
            items: string
        inputBinding:
            position: 14
            itemSeparator: ' , '
            shellQuote: False
            prefix: '-a'
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
            prefix: '-h'
    chim_nonchim_scoredrop_min:
        type: int
        default: 10
        inputBinding:
            position: 16
            prefix: '-d'
    peoverlap_nbases_min:
        type: int
        default: 12
        inputBinding:
            position: 17
            prefix: '-o'
    peoverlap_mmp:
        type: float
        default: 0.1
        inputBinding:
            position: 18
            prefix: '-p'
    chimout_junction_format:
        type: int
        default: 1
        inputBinding:
            position: 19
            prefix: '-v'
    star_genome_dir:
        type: Directory
        inputBinding:
            position: 20
            prefix: '-g'
        doc: 'specifies path to the directory where the genome indices are stored'
    twopass_mode:
        type: string
        default: "Basic"
        inputBinding:
            position: 21
            prefix: '-w'
    gtf_file:
        type: File?
        inputBinding:
            position: 22
            prefix: '-k'
    paired_end:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
        inputBinding:
            position: 26
            prefix: '-e'
    fastqs:
        type:
            type: array
            items:
                type: array
                items: File
        inputBinding:
            position: 27
            prefix: '-f'

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
