#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Helper script to pull the k-mer size used in generating a kallisto index from the index file"
baseCommand: ["perl", "extractor.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "ubuntu:bionic"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'extractor.pl'
        entry: |
            my $filename = $ARGV[0];
            open(my $fh, $filename);
            binmode $fh;
            read($fh, my $version, 8);
            die "unexpected index version!" unless unpack('Q', $version) == 10;
            read($fh, my $kmer_size, 4);
            print unpack('L', $kmer_size);
stdout: kmer_temp.txt
inputs:
    kallisto_index:
        type: File
        inputBinding:
            position: 1
outputs:
    kmer_length:
        type: int
        outputBinding:
            glob: "kmer_temp.txt"
            loadContents: true
            outputEval: $(parseInt(self[0].contents))
