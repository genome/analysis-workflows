#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: "STAR: align reads to transcriptome"
baseCommand: ["/usr/local/bin/STAR"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 10
    - class: DockerRequirement
      dockerPull: "sridnona/mgirnaseq:v3"
arguments: [
    "--runThreadN", $(runtime.cores)

]
inputs:
    outSAMtype:
        type: string[]
        default: ["BAM", "Unsorted"]
        inputBinding:
            position: 3
            prefix: '--outSAMtype'
    runMode:
        type: string
        default: "alignReads"
        inputBinding:
            position: 2
            prefix: "--runMode"
        doc: |
          string: type of the run
          alignReads - map reads
          genomeGenerate - generate genome files
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
    outReadsUnmapped:
        type: string
        default: "None"
        inputBinding:
            position: 6
            prefix: '--outReadsUnmapped'
    chimSegmentMin:
        type: int
        default: "12"
        inputBinding:
            position: 7
            prefix: '--chimSegmentMin'
    chimJunctionOverhangMin:
        type: int
        default: "12"
        inputBinding:
            position: 8
            prefix: '--chimJunctionOverhangMin'
    alignSJDBoverhangMin:
        type: int
        default: 10
        inputBinding:
            position: 9
            prefix: '--alignSJDBoverhangMin'
    alignMatesGapMax:
        type: int?
        default: 100000
        inputBinding:
            position: 10
            prefix: '--alignMatesGapMax'
    alignIntronMax:
        type: int
        default: 100000
        inputBinding:
            position: 11
            prefix: '--alignIntronMax'
    chimSegmentReadGapMax:
        type: int
        default: 3
        inputBinding:
            position: 12
            prefix: '--chimSegmentReadGapMax'
    alignSJstitchMismatchNmax:
        type: int[]
        default: [5, -1, 5, 5]
        itemSeparator: ' '
        inputBinding:
            position: 13 
            prefix: '--alignSJstitchMismatchNmax'
    outSAMstrandField:
        type: string
        default: "intronMotif"
        inputBinding:
            position: 14
            prefix: '--outSAMstrandField'
    outSAMunmapped:
        type: string
        default: Within
        inputBinding:
            position: 15
            prefix: '--outSAMunmapped'
    outSAMattrRGline:
        type:
            type: array
            items: string
        inputBinding:
            position: 16
            itemSeparator: ' , ' 
            shellQuote: False
            prefix: '--outSAMattrRGline'
    chimMultimapNmax:
        type: int
        default: 10
        inputBinding:
            position: 17
            prefix: '--chimMultimapNmax'
    chimNonchimScoreDropMin:
        type: int
        default: 10
        inputBinding:
            position: 18
            prefix: '--chimNonchimScoreDropMin'
    peOverlapNbasesMin:
        type: int
        default: 12
        inputBinding:
            position: 19
            prefix: '--peOverlapNbasesMin'
    peOverlapMMp:
        type: float
        default: 0.1
        inputBinding:
            position: 20
            prefix: '--peOverlapMMp'
    chimOutJunctionFormat:
        type: int
        default: 1
        inputBinding:
            position: 21
            prefix: '--chimOutJunctionFormat'
    stargenomeDir:
        type: Directory
        inputBinding:
            position: 22
            prefix: '--genomeDir'
    twopassMode:
        type: string
        default: "Basic"
        inputBinding:
            position: 23
            prefix: '--twopassMode'
    gtf_file:
        type: File?
        inputBinding:
            position: 24
            prefix: '--sjdbGTFfile'
    outFileNamePrefix:
        type: string
        default: "STAR_"
        inputBinding:
            position: 25
            prefix: '--outFileNamePrefix'
    readFilesCommand:
        default: "cat"
        type: string?
        inputBinding:
            position: 26
            prefix: '--readFilesCommand'
    outSAMattributes:
        type: string[]
        default: [NH, HI, AS, NM, MD]
        itemSeparator: ' '
        inputBinding:
            position: 27
            prefix: '--outSAMattributes'

outputs:
    aligned_bam:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)Aligned.out.bam"
    log_final:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)Log.final.out"
    log:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)Log.out"
    log_progress:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)Log.progress.out"
    splicejunctionout:
        type: File
        outputBinding:
            glob: "$(inputs.outFileNamePrefix)SJ.out.tab"
    chimjunc:
        type: File
        outputBinding:
            glob: "$(inputs.outFileNamePrefix)Chimeric.out.junction"
