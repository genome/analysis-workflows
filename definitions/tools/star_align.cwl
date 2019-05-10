#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0
label: "STAR: align reads to transcriptome"
baseCommand: STAR
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 16
    - class: DockerRequirement
      dockerPull: "sridnona/mgirnaseq:v1"
arguments: [
    "--runThreadN", $(runtime.cores)

]
inputs:
    outSAMtype:
        type: string[]
        default: ["BAM", "SortedByCoordinate"]
        inputBinding:
            position: 3
            prefix: '--outSAMtype'
    runMode:
        type: string
        default: alignReads
        inputBinding:
            position: 2
            prefix: "--runMode"
        doc: |
          string: type of the run
          alignReads             ... map reads
          genomeGenerate         ... generate genome files
          inputAlignmentsFromBAM ... input alignments from BAM. Presently only works with alignReads.
    fastq1:
        type: File
        inputBinding:
            position: 4
            prefix: '--readFilesIn'
    fastq2:
        type: File
        inputBinding:
            position: 5
            prefix: ''
    stargenomeDir:
        type: Directory
        inputBinding:
            position: 6
            prefix: '--genomeDir'
        doc: |
          string: paths to the directory where genome files are stored, created using genomeGenerate
    outSAMunmapped:
        type: int?
        inputBinding:
            position: 7
            prefix: '--outSAMunmapped'
        doc: |
          string: output of unmapped reads in the SAM format
          None   ... no output
          Within ... output unmapped reads within the main SAM file (i.e. Aligned.out.sam)
    alignIntronMin:
        type: int?
        inputBinding:
            position: 8
            prefix: '--alignIntronMin'
        doc: 'minimum intron size: genomic gap is considered intron if its length>=alignIntronMin,otherwise it is considered Deletion
        '
    alignIntronMax:
        type: int?
        inputBinding:
            position: 9
            prefix: '--alignIntronMax'
        doc: 'maximum intron size: genomic gap is considered intron if its length>=alignIntronMin,otherwise it is considered Deletion
        '
    alignMatesGapMax:
        type: int?
        inputBinding:
            position: 10
            prefix: '--alignMatesGapMax'
        doc: 'maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins
        '
    outSAMstrandField:
        type: string
        default: intronMotif
        inputBinding:
            position: 11
            prefix: '--outSAMstrandField'
        doc: |
          None
          string: Cufflinks-like strand field flag
          None        ... not used
          intronMotif ... strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.
          NOTE you need this option to be set yto intronmotif if for stringtie quantification
    twopassMode:
        type: string?
        inputBinding:
            position: 12
            prefix: '--twopassMode'
        doc: |
           None
           string: 2-pass mapping mode.
           None        ... 1-pass mapping
           Basic ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly
    outFilterMultimapNmax:
        type: int?
        inputBinding:
            position: 13
            prefix: '--outFilterMultimapNmax'
        doc: 'int: read alignments will be output only if the read maps fewer than this value,otherwise no alignments will be output
        '
    sjdbOverhang:
        type: string?
        inputBinding:
            position: 14
            prefix: '--sjdbOverhang'
        doc: 'int>0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)
        '
    gtf_file:
        type: File?
        inputBinding:
            position: 15
            prefix: '--sjdbGTFfile'
        doc: |
          string: path to the GTF file with annotations
    outFileNamePrefix:
        type: string
        inputBinding:
            position: 16
            prefix: '--outFileNamePrefix'
        doc: ' output files name prefix (including full or relative path). Can
          only be defined on the command line 
          '
    readFilesCommand:
        type: string?
        inputBinding:
            position: 17
            prefix: '--readFilesCommand'
        doc: ' command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout, zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc 
        '
    outSAMattributes:
        type: string[]
        default: ["NH", "HI", "AS", "NM", "MD"]
        inputBinding:
            position: 18
            prefix: '--outSAMattributes'
        doc: |
          Standard
          string: a string of desired SAM attributes, in the order desired for the output SAM
          NH HI AS nM NM MD jM jI XS ... any combination in any order
          Standard   ... NH HI AS nM
          All        ... NH HI AS nM NM MD jM jI
          None ... no attributes
    outSAMattrRGline:
        type:
            type: array
            items: string
        inputBinding:
            position: 19
            prefix: '--outSAMattrRGline'
        doc: |
           -
           string(s): SAM/BAM read group line. The first word contains the read group identifier and must start with "ID:", e.g. --outSAMattrRGline ID:xxx CN:yy "DS:z z z".
           xxx will be added as RG tag to each output alignment. Any spaces in the tag values have to be double quoted.
           Comma separated RG lines correspons to different (comma separated) input files in --readFilesIn. Commas have to be surrounded by spaces, e.g.
           --outSAMattrRGline ID:xxx , ID:zzz "DS:z z" , ID:yyy DS:yyyy


outputs:
    aligned_bam:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)_Aligned.sortedByCoord.out.bam"
    log_final:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)_Log.final.out"
    log:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)_Log.out"
    log_progress:
        type: File
        outputBinding:
          glob: "$(inputs.outFileNamePrefix)_Log.progress.out"
    splicejunctionout:
        type: File
        outputBinding:
            glob: "$(inputs.outFileNamePrefix)_SJ.out.tab"

