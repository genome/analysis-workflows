#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 somatic"
baseCommand: "varscan_helper.sh"
requirements:
    - class: ResourceRequirement
      ramMin: 12000
      coresMin: 2
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.3.1"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'varscan_helper.sh'
        entry: |
            #!/bin/bash

            set -o errexit
            set -o nounset

            if [ $# -lt 7 ]
            then
                echo "Usage: $0 [TUMOR_BAM] [NORMAL_BAM] [REFERENCE] [STRAND_FILTER] [MIN_COVERAGE] [MIN_VAR_FREQ] [P_VALUE] [roi_bed?]"
                exit 1
            fi

            TUMOR_BAM="$1"
            NORMAL_BAM="$2"
            REFERENCE="$3"
            STRAND_FILTER="$4"
            MIN_COVERAGE="$5"
            MIN_VAR_FREQ="$6"
            P_VALUE="$7"
            OUTPUT="${HOME}/output"

            if [ -z ${8+x} ]
            then
                #run without ROI
                java -jar /opt/varscan/VarScan.jar somatic \
                    <(/opt/samtools/bin/samtools mpileup --no-baq -f "$REFERENCE" "$NORMAL_BAM" "$TUMOR_BAM") \
                    $OUTPUT \
                    --strand-filter $STRAND_FILTER \
                    --min-coverage $MIN_COVERAGE \
                    --min-var-freq $MIN_VAR_FREQ \
                    --p-value $P_VALUE \
                    --mpileup 1 \
                    --output-vcf
            else
                ROI_BED="$8"
                java -jar /opt/varscan/VarScan.jar somatic \
                    <(/opt/samtools/bin/samtools mpileup --no-baq -l "$ROI_BED" -f "$REFERENCE" "$NORMAL_BAM" "$TUMOR_BAM") \
                    $OUTPUT \
                    --strand-filter $STRAND_FILTER \
                    --min-coverage $MIN_COVERAGE \
                    --min-var-freq $MIN_VAR_FREQ \
                    --p-value $P_VALUE \
                    --mpileup 1 \
                    --output-vcf
            fi

inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [^.bai]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            position: 3
    strand_filter:
        type: int?
        default: 0
        inputBinding:
            position: 4
    min_coverage:
        type: int?
        default: 8
        inputBinding:
            position: 5
    min_var_freq:
        type: float?
        default: 0.1
        inputBinding:
            position: 6
    p_value:
        type: float?
        default: 0.99
        inputBinding:
            position: 7
    roi_bed:
        type: File?
        inputBinding:
            position: 8
outputs:
    snvs:
        type: File
        outputBinding:
            glob: "output.snp.vcf"
    indels:
        type: File
        outputBinding:
            glob: "output.indel.vcf"
