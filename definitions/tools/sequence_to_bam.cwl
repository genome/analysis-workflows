#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "resolve sequence type to a bam"
baseCommand: ["/bin/bash", "sequence_to_bam_helper.sh"]
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "broadinstitute/picard:2.23.6"
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'sequence_to_bam_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit
            set -o nounset

            if [[ "$1" == 'bam' ]]; then
                cp "$2" sequence.bam
            else
                /usr/bin/java -Xmx4g -jar /usr/picard/picard.jar FastqToSam --OUTPUT sequence.bam "$@"
            fi
arguments:
    - { valueFrom: "$(inputs.sequence.sequence.hasOwnProperty('bam')? inputs.sequence.sequence.bam : null)", prefix: 'bam', position: -1 }
    - { valueFrom: "$(inputs.sequence.sequence.hasOwnProperty('fastq1')? inputs.sequence.sequence.fastq1 : null)", prefix: '--FASTQ' }
    - { valueFrom: "$(inputs.sequence.sequence.hasOwnProperty('fastq2')? inputs.sequence.sequence.fastq2 : null)", prefix: '--FASTQ2' }
    - valueFrom: |
          ${
              var x = inputs.sequence.readgroup.split("\t").find(function(tag){ return tag.startsWith("SM:")});
              if(x) {
                  return x.substr(3)
              } else {
                  return null;
              }
          }
      prefix: '--SAMPLE_NAME'
    - valueFrom: |
          ${
              var x = inputs.sequence.readgroup.split("\t").find(function(tag){ return tag.startsWith("LB:")});
              if(x) {
                  return x.substr(3)
              } else {
                  return null;
              }
          }
      prefix: '--LIBRARY_NAME'
    - valueFrom: |
          ${
              var x = inputs.sequence.readgroup.split("\t").find(function(tag){ return tag.startsWith("PU:")});
              if(x) {
                  return x.substr(3)
              } else {
                  return null;
              }
          }
      prefix: '--PLATFORM_UNIT'
    - valueFrom: |
          ${
              var x = inputs.sequence.readgroup.split("\t").find(function(tag){ return tag.startsWith("PL:")});
              if(x) {
                  return x.substr(3)
              } else {
                  return null;
              }
          }
      prefix: '--PLATFORM'
    - valueFrom: |
          ${
              var x = inputs.sequence.readgroup.split("\t").find(function(tag){ return tag.startsWith("ID:")});
              if(x) {
                  return x.substr(3)
              } else {
                  return null;
              }
          }
      prefix: '--READ_GROUP_NAME'
    - valueFrom: |
          ${
              var x = inputs.sequence.readgroup.split("\t").find(function(tag){ return tag.startsWith("CN:")});
              if(x) {
                  return x.substr(3)
              } else {
                  return null;
              }
          }
      prefix: '--SEQUENCING_CENTER'
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data
    output_bam_basename:
        type: string
        default: "sequence"
outputs:
    bam:
        type: File
        outputBinding:
            glob: sequence.bam
