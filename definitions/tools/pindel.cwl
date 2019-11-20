#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "pindel v0.2.5b8"
arguments: ["/usr/bin/perl", "pindel_helper.pl"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      tmpdirMin: 100000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.4.2"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'pindel_helper.pl'
        entry: |
            #!/usr/bin/perl

            use strict;
            use warnings;

            use IO::File;

            unless (@ARGV > 5) {
                die "Usage: $0 normal.bam tumor.bam insert_size normal_sample_name tumor_sample_name <args>";
            }

            my ($normal_bam, $tumor_bam, $insert_size, $normal_name, $tumor_name, @args) = @ARGV;

            my $fh = IO::File->new("> pindel.config");

            $fh->say(join("\t", $normal_bam, $insert_size, $normal_name));
            $fh->say(join("\t", $tumor_bam, $insert_size, $tumor_name));
            $fh->close;

            exit system(qw(/usr/bin/pindel -i pindel.config -w 30 -T 4 -o all), @args);

inputs:
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
        inputBinding:
            position: 1
    normal_bam:
        type: File
        secondaryFiles: ["^.bai"]
        inputBinding:
            position: 2
    insert_size:
        type: int
        default: 400
        inputBinding:
            position: 3
    tumor_sample_name:
        type: string
        inputBinding:
            position: 4
    normal_sample_name:
        type: string
        inputBinding:
            position: 5
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-f"
            position: 6
    chromosome:
        type: string?
        inputBinding:
            prefix: "-c"
            position: 7
    region_file:
        type: File?
        inputBinding:
            prefix: "-j"
            position: 8
outputs:
    deletions:
        type: File
        outputBinding:
            glob: "all_D"
    insertions:
        type: File
        outputBinding:
            glob: "all_SI"
    tandems:
        type: File
        outputBinding:
            glob: "all_TD"
    long_insertions:
        type: File
        outputBinding:
            glob: "all_LI"
    inversions:
        type: File
        outputBinding:
            glob: "all_INV"
