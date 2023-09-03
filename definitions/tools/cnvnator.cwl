#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run CNVnator to calculate copy number variations in WGS samples"

baseCommand: ["/bin/bash", "run_cnvnator.sh"]

requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cnvnator-cwl:0.4"
    - class: ResourceRequirement
      ramMin: 20000
      coresMin: 1
      tmpdirMin: 10000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "run_cnvnator.sh"
        entry: |
          #!/bin/bash

          #set up the environment
          source /opt/root/bin/thisroot.sh

          set -eou pipefail

          # set vars
          BAM="$1"
          BIN_SIZE="$2"
          CHROMOSOMES="$3"
          REFERENCE="$4"
          SAMPLE="$5"

          # create directory to store fasta files. CNVnator wants chromosomes in individual files
          # while with later versions of CNVnator these steps will accept a fasta.gz file the cnvnator2VCF.pl will not
          mkdir FASTA_CHRS
          awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > "FASTA_CHRS/"CHROM".fa" }' "$REFERENCE"

          # extract read mapping from input bam(single sample)
          cnvnator -root "$SAMPLE.root" -tree "$BAM" -chrom $CHROMOSOMES
          # generate read depth histogram
          cnvnator -root "$SAMPLE.root" -his "$BIN_SIZE" -d FASTA_CHRS/ -chrom $CHROMOSOMES
          # calculate statistics
          cnvnator -root "$SAMPLE.root" -stat "$BIN_SIZE" -chrom $CHROMOSOMES
          # read depth signal partitioning
          cnvnator -root "$SAMPLE.root" -partition "$BIN_SIZE" -chrom $CHROMOSOMES
          # cnv calling
          cnvnator -root "$SAMPLE.root" -call "$BIN_SIZE" -chrom $CHROMOSOMES > "$SAMPLE.cnvnator.cn"

          # convert to vcf
          cnvnator2VCF.pl -reference "$REFERENCE" "$SAMPLE.cnvnator.cn" FASTA_CHRS/ >  "$SAMPLE.cnvnator.vcf"
          exit 0
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
        doc: "BAM or CRAM input"
    bin_size:
        type: int?
        default: 100
        inputBinding:
            position: 2
        doc: "Bin size used to calculate coverage"
    chromosomes:
        type: string[]?
        default: ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        inputBinding:
           position: 3
           itemSeparator: " "
        doc: "List of chromosomes to run CNVnator on"
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai]
        inputBinding:
            position: 4
        doc: "Reference used to generate the alignments"
    sample_name:
        type: string
        inputBinding:
            position: 5
        doc: "Sample name used for output file naming"

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).cnvnator.vcf"
    root_file:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).root"
    cn_file:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).cnvnator.cn"
