#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', '-M5.10.0', '-ne', 'if(substr($_,0,1) ne q{@}) { chomp; my @c = split "\t"; say(join("\t", $c[0], $c[1]-1, $c[2])); }']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdin: $(inputs.interval_list.path)
stdout: "interval_list.bed"
inputs:
    interval_list:
        type: File
outputs:
    interval_bed:
        type: stdout

