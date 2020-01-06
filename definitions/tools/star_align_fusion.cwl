#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: "STAR: align reads to transcriptome"
baseCommand: ["/bin/bash","star.sh"]
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

            while getopts "c:f:s:j:b:q:i:r:m:t:u:a:h:d:o:p:v:g:w:k:e:" opt; do
                case "$opt" in
                    c) cores="$OPTARG";;
                    f) fastqs="$OPTARG";;
                    s) chimSegmentMin="$OPTARG";;
                    j) chimJunctionOverhangMin="$OPTARG";;
                    b) alignSJDBoverhangMin="$OPTARG";;
                    q) alignMatesGapMax="$OPTARG";;
                    i) alignIntronMax="$OPTARG";;
                    r) chimSegmentReadGapMax="$OPTARG";;
                    m) alignSJstitchMismatchNmax="$OPTARG";;
                    u) outSAMunmapped="$OPTARG";;
                    a) outSAMattrRGline="$OPTARG";;
                    h) chimMultimapNmax="$OPTARG";;
                    d) chimNonchimScoreDropMin="$OPTARG";;
                    o) peOverlapNbasesMin="$OPTARG";;
                    p) peOverlapMMp="$OPTARG";;
                    v) chimOutJunctionFormat="$OPTARG";;
                    g) genomeDir="$OPTARG";;
                    w) twopassMode="$OPTARG";;
                    e) paired="$OPTARG";;
                    k) sjdbGTFfile="$OPTARG";;
                esac
            done
            echo $cores
            echo $fastqs
            echo $chimSegmentMin
            echo $chimJunctionOverhangMin
            echo $alignSJDBoverhangMin
            echo $alignMatesGapMax
            echo $alignIntronMax
            echo $chimSegmentReadGapMax
            echo $alignSJstitchMismatchNmax
            echo $outSAMunmapped
            echo $outSAMattrRGline
            echo $chimMultimapNmax
            echo $chimNonchimScoreDropMin
            echo $peOverlapNbasesMin
            echo $peOverlapMMp
            echo $chimOutJunctionFormat
            echo $genomeDir
            echo $twopassMode
            echo $sjdbGTFfile
            echo $paired

            fqfinal=""
            if [[ "$paired" == "false" ]];then
                #run star with all the fastqs in single-end mode
                fqfinal=$fastqs
            else
                #split interleaved fastqs into fq1/2 arguments
                #can't use arrays so this gets a bit unweildy
                fq1=""
                fq2=""
                i=0
                oldifs=$IFS
                IFS=$'\n' # can't use standard separator, because spaces could have files in the name
                # this gets a bit unweildy because array syntax 
                # chokes CWL (can't use dollar sign curly braces)
                for fq in `echo $fastqs | tr "," "\n"`;do
                    echo "test $fq"
                    mod=`echo $i % 2`
                    if [[ $mod -eq 0 ]];then
                        count=`printf "%s" "$fq1" | wc -c`
                        if [[ $count -gt 0 ]];then
                            fq1="$fq1,$fq";
                        else
                            fq1="$fq";
                        fi
                    else
                        count=`printf "%s" "$fq2" | wc -c`
                        if [[ $count -gt 0 ]];then
                            fq2="$fq2,$fq";
                        else
                            fq2="$fq";
                        fi
                    fi
                    let "i=$i+1"
                 done
                IFS=$oldifs #clean up after ourselves (though I don't think it's really needed)                
                fqfinal="$fq1 $fq2"
            fi

            echo "/usr/local/bin/STAR --runMode alignReads --outSAMtype BAM Unsorted --outReadsUnmapped None --outFileNamePrefix STAR_ --readFilesCommand cat --outSAMattributes NH HI AS NM MD --runThreadN \"$cores\" --readFilesIn \"$fqfinal\" --chimSegmentMin \"$chimSegmentMin\" --chimJunctionOverhangMin \"$chimJunctionOverhangMin\" --alignSJDBoverhangMin \"$alignSJDBoverhangMin\" --alignMatesGapMax \"$alignMatesGapMax\" --alignIntronMax \"$alignIntronMax\" --chimSegmentReadGapMax \"$chimSegmentReadGapMax\" --alignSJstitchMismatchNmax \"$alignSJstitchMismatchNmax\" --outSAMstrandField \"intronMotif\" --outSAMunmapped \"$outSAMunmapped\" --outSAMattrRGline \"$outSAMattrRGline\" --chimMultimapNmax \"$chimMultimapNmax\" --chimNonchimScoreDropMin \"$chimNonchimScoreDropMin\" --peOverlapNbasesMin \"$peOverlapNbasesMin\" --peOverlapMMp \"$peOverlapMMp\" --chimOutJunctionFormat \"$chimOutJunctionFormat\" --genomeDir \"$genomeDir\" --twopassMode \"$twopassMode\" --sjdbGTFfile \"$sjdbGTFfile\""

            /usr/local/bin/STAR --runMode alignReads --outSAMtype BAM Unsorted --outReadsUnmapped None --outFileNamePrefix STAR_ --readFilesCommand cat --outSAMattributes NH HI AS NM MD --runThreadN $cores --readFilesIn $fqfinal --chimSegmentMin $chimSegmentMin --chimJunctionOverhangMin $chimJunctionOverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --alignMatesGapMax $alignMatesGapMax --alignIntronMax $alignIntronMax --chimSegmentReadGapMax $chimSegmentReadGapMax --alignSJstitchMismatchNmax $alignSJstitchMismatchNmax --outSAMstrandField intronMotif --outSAMunmapped $outSAMunmapped --outSAMattrRGline $outSAMattrRGline --chimMultimapNmax $chimMultimapNmax --chimNonchimScoreDropMin $chimNonchimScoreDropMin --peOverlapNbasesMin $peOverlapNbasesMin --peOverlapMMp $peOverlapMMp --chimOutJunctionFormat $chimOutJunctionFormat --genomeDir $genomeDir --twopassMode $twopassMode --sjdbGTFfile $sjdbGTFfile

#touch STAR_Aligned.out.bam STAR_Log.final.out STAR_Log.out STAR_Log.progress.out STAR_SJ.out.tab STAR_Chimeric.out.junction

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
            shellQuote: True
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
            shellQuote: True
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
        type: File[]
        inputBinding:
            position: 27
            prefix: "-f"
            itemSeparator: ","
outputs:
    aligned_bam:
        type: File
        outputBinding:
          glob: "STAR_Aligned.out.bam"
    log_final:
        type: File
        outputBinding:
          glob: "STAR_Log.final.out"
    log:
        type: File
        outputBinding:
          glob: "STAR_Log.out"
    log_progress:
        type: File
        outputBinding:
          glob: "STAR_Log.progress.out"
    splice_junction_out:
        type: File
        outputBinding:
            glob: "STAR_SJ.out.tab"
    chim_junc:
        type: File
        outputBinding:
            glob: "STAR_Chimeric.out.junction"
