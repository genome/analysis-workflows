#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "merges nearby DEL/DUP records within a certain window distance"

baseCommand: ["python3", "merge.py"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "griffithlab/vatools:4.1.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: "merge.py"
        entry: |
          import argparse
          import vcfpy
          from collections import OrderedDict

          parser = argparse.ArgumentParser()
          parser.add_argument('--input', '-i', dest="input", help='input vcf file', required=True, action="store")
          parser.add_argument('--output', '-o', dest="output", help='output vcf file', required=False, default="out.vcf", action="store")
          parser.add_argument('--window', '-w', dest="window", help='max merge window size', required=False, default=1000, type=int, action="store")

          args = parser.parse_args()
          in_vcf_name  = args.input
          out_vcf_name = args.output
          window_size = args.window

          reader = vcfpy.Reader.from_path(in_vcf_name)
          new_header = reader.header
          new_header.add_filter_line(vcfpy.OrderedDict([('ID', 'MERGED_CALL'), ('Description', 'Record merged from 2 or more individual records')]))

          writer = vcfpy.Writer.from_path(out_vcf_name, new_header)
          new_record_count = 0
          merge_records = []
          for record in reader:
              if((len(merge_records) == 0) or (merge_records[-1].CHROM != record.CHROM) or (merge_records[-1].INFO['SVTYPE'] != record.INFO['SVTYPE']) or (abs(merge_records[-1].INFO['END'] - record.POS) > window_size)):
                  if(len(merge_records) > 1):
                      new_record_count = new_record_count + 1
                      new_record_chr = merge_records[0].CHROM
                      new_record_start = merge_records[0].POS
                      new_record_end = merge_records[-1].INFO['END']
                      new_record_type = merge_records[0].INFO['SVTYPE']
                      new_record_svlen = new_record_end - new_record_start

                      info = OrderedDict({"SVTYPE": new_record_type, "END": new_record_end, "SVLEN": new_record_svlen})
                      alt = vcfpy.SymbolicAllele(new_record_type)
                      sample_calls = []
                      for sample in merge_records[0].calls:
                          gt = OrderedDict({"GT": "/".join(map(str, sample.gt_alleles)).replace("None",".")})
                          name = sample.sample
                          sample_calls.append(vcfpy.Call(name, gt))

                      new_record = vcfpy.Record(new_record_chr, new_record_start, [], "N", [alt], ".", ["MERGED_CALL"], info, ["GT"], sample_calls)
                      writer.write_record(new_record)
                      merge_records = [record]
                  else:
                      merge_records = [record]
                  writer.write_record(record)
                  next

              dist = abs(merge_records[-1].INFO['END'] - record.POS)
              if(dist < window_size):
                  merge_records.append(record)
          print(f"Found {new_record_count} records that can be merged based on the input {window_size} distance")

inputs:
    input_vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    output_vcf_name:
        type: string?
        default: "record_merged.vcf"
        inputBinding:
            prefix: "-o"
            position: 2
    distance:
        type: int?
        default: 1000
        inputBinding:
            prefix: "-w"
            position: 3

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "$(inputs.output_vcf_name)"
