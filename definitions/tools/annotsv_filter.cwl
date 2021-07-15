#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "filter AnnotSV tsv output generated from a bcftools merged vcf"

baseCommand: ["python3", "AnnotSV_filter.py"]

requirements:
    - class: DockerRequirement
      dockerPull: "python:3"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "AnnotSV_filter.py"
        entry: |
            import argparse
            import csv
            import sys

            parser = argparse.ArgumentParser()
            parser.add_argument('--input', '-i', dest="input", help='input AnnotSV tsv file', required=True, action="store")
            parser.add_argument('--output', '-o', dest="output", help='output tsv file name', required=True, action="store")
            parser.add_argument('--filtering_frequency', dest="filtering_frequency", help="frequency to filter with", action="store", type=float, default="0.05")
            parser.add_argument('--no-CDS', dest="CDS", help="Do not require a positive CoDing Sequence overlap", action="store_true")
            parser.add_argument('--ignore-pass-filter', dest="filter", help="Do not require calls to have a PASS filter", action="store_true")

            args = parser.parse_args()
            input_file_name  = args.input
            output_file_name = args.output
            filtering_frequency = args.filtering_frequency
            no_cds = args.CDS
            ignore_pass_filter = args.filter

            with open(input_file_name, 'r') as file_in, open(output_file_name, 'w') as file_out:
                file_in = csv.DictReader(file_in, delimiter='\t')
                file_out = csv.DictWriter(file_out, fieldnames=file_in.fieldnames, delimiter='\t')
                file_out.writeheader()
                total_sv_count = 0
                pass_sv_count = 0
                for row in file_in:
                    total_sv_count += 1
                    if(row['AnnotSV type'] == 'split' \
                        and (row['FILTER'] == 'PASS' or ignore_pass_filter) \
                        and (int(row['CDS length']) > 0 or no_cds) \
                        and float(row['IMH_AF']) < filtering_frequency
                        and float(row['1000g_max_AF']) < filtering_frequency
                        and not(float(row['DGV_LOSS_Frequency']) > filtering_frequency and 'DEL' in row['SV type']) 
                        and not(float(row['DGV_GAIN_Frequency']) < filtering_frequency and ('DUP' in row['SV type'] or 'INS' in row['SV type']))
                        and not(('Manta' in row['ID'] and 'IMPRECISE' in row['INFO']) or (row['QUAL'] != '.' and 'IMPRECISE' in row['INFO'])) ):
                        file_out.writerow(row)
                        pass_sv_count += 1
                print("total sv count:",total_sv_count)
                print("total sv passed count:",pass_sv_count)

inputs:
    no_CDS:
        type: boolean?
        inputBinding:
            position: 1
            prefix: "--no-CDS"
    annotsv_tsv:
        type: File
        inputBinding:
            position: 2
            prefix: "--input"
    filtering_frequency:
        type: double?
        inputBinding:
            position: 3
            prefix: "--filtering_frequency"
    ignore_pass_filter:
        type: boolean?
        inputBinding:
            position: 4
            prefix: "--ignore-pass-filter"
    output_tsv_name:
        type: string?
        default: "filtered-bcftools-merged-AnnotSV.tsv"
        inputBinding:
            position: 5
            prefix: "--output"

outputs:
    filtered_tsv:
        type: File
        outputBinding:
            glob: $(inputs.output_tsv_name)
