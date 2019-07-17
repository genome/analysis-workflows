#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Helper script to pull the k-mer size used in generating a kallisto index from the index file"
baseCommand: ["g++", "-o", "extractor", "extractor.cpp"]
arguments: [{ shellQuote: false, valueFrom: "&&" }, "./extractor"]
requirements:
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 16000
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq:latest"
    - class: ShellCommandRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'extractor.cpp'
        entry: |
           //based on https://github.com/pachterlab/kallisto/blob/master/src/KmerIndex.cpp#L809
           //accessed 7/11/2019
           
           #include <string>
           #include <iostream>
           #include <fstream>
           
           //single argument is the index file from which k-mer size will be extracted
           int main(int argc, char *argv[]) {
           
             std::string index(argv[1]);
             const char *index_in = index.c_str();
             std::ifstream in(index_in, std::ios::in | std::ios::binary);
           
             //ignore first part of header
             in.ignore(sizeof(size_t), EOF);

             // read k
             int k;
             in.read((char *)&k, sizeof(k));

             std::ofstream outfile;
             outfile.open("helper.txt");
             outfile << k;
             outfile.close();
             return 0;
           
           }

inputs:
    kallisto_index:
        type: File
        inputBinding:
            position: 1
outputs:
    kmer_length:
        type: int
        outputBinding:
            glob: "helper.txt"
            loadContents: true
            outputEval: $(parseInt(self[0].contents))
