#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: "Script to create FDA-requested summary tables"
requirements:
    - class: DockerRequirement
      dockerPull: "python:3.7.4-slim-buster"
    - class: ResourceRequirement
      ramMin: 8000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'generate_tables.py'
        entry: |
            import argparse, zipfile, os

            parser = argparse.ArgumentParser()

            #data files to be parsed
            parser.add_argument('--table', choices=['table1', 'table2', 'table3'], required=True)
            parser.add_argument('--md5sum_file', required=True)
            parser.add_argument('--fastqc_zips', nargs='*', required=True)
            parser.add_argument('--aligned_metrics')
            parser.add_argument('--unaligned_metrics', nargs='*')
            parser.add_argument('--alignment_summary_metrics')
            parser.add_argument('--duplication_metrics')
            parser.add_argument('--insert_size_metrics')
            parser.add_argument('--hs_metrics')
            parser.add_argument('--rna_metrics')
            parser.add_argument('--flagstat')
            parser.add_argument('--unaligned_rna_table')

            #inputs to immuno pipeline
            parser.add_argument('--sequencing_platform', default='NOT_PROVIDED')
            parser.add_argument('--sequencing_instrument', default='NOT_PROVIDED')
            parser.add_argument('--sequencing_kit', default='NOT_PROVIDED')
            parser.add_argument('--sequencing_type', default='NOT_PROVIDED')
            parser.add_argument('--single_or_paired_end', default='NOT_PROVIDED')
            parser.add_argument('--spike_in_error_rate', default='NOT_PROVIDED')
            parser.add_argument('--source', default='NOT_PROVIDED')
            parser.add_argument('--total_DNA', default='NOT_PROVIDED')
            parser.add_argument('--reference_genome', default='NOT_PROVIDED')
            parser.add_argument('--total_RNA', default='NOT_PROVIDED')
            parser.add_argument('--RIN_score', default='NOT_PROVIDED')
            parser.add_argument('--freq_normalization_method', default='NOT_PROVIDED')
            parser.add_argument('--annotation_file', default='NOT_PROVIDED')
            parser.add_argument('--sample_name', default='NOT_PROVIDED')

            parser.add_argument('--table_file_name', required=True)

            args = parser.parse_args()

            # map table row names to values pulled directly from command line inputs
            args_to_table_fields = {
                "Sequencing Platform": args.sequencing_platform,
                "Sequencing Instrument": args.sequencing_instrument,
                "Sequencing Kit": args.sequencing_kit,
                "Sequencing Type": args.sequencing_type,
                "Single or Paired End": args.single_or_paired_end,
                "Spike in Error Rate (%)": args.spike_in_error_rate,
                "Source": args.source,
                "Total DNA(ng)": args.total_DNA,
                "Reference Genome": args.reference_genome,
                "Total RNA (ng)": args.total_RNA,
                "RIN score": args.RIN_score,
                "Frequency Normalization Method": args.freq_normalization_method,
                "Annotation File": args.annotation_file,
                "Sample Name": args.sample_name,
                "FastQC Plots Provided?": 'Yes'
            }

            # the following 2 methods produce nested dictionaries, with the top level keys
            # corresponding to file names and their sub-dictionaries containing table row
            # names and values

            def parse_md5sums(md5_file):
                with open(md5_file) as f:
                    parsed_data = [line.split() for line in f.read().splitlines()]
                extracted_data = {}
                for md5sum, filepath in parsed_data:
                    filename = os.path.basename(filepath)
                    extracted_data[filename] = {"MD5 File Checksum": md5sum}
                    extracted_data[filename]['File'] = filename
                return extracted_data

            def parse_fastqc_zips(zip_list, table_num):

                #define mappings from field name in fastqc data to row label in table
                fieldname_to_t1_label = {
                    "Total Sequences": "Total Number of Reads",
                    "%GC": "GC Content (%)"
                }

                fieldname_to_t2_label = {
                    "Sequence length": "Read Length"
                }

                fieldname_to_t3_label = {
                    "Sequence length": "Read Length (nt)", 
                }

                #helper to choose proper mapping dict for requested table
                table_to_map_dict = {
                    "table1": fieldname_to_t1_label,
                    "table2": fieldname_to_t2_label,
                    "table3": fieldname_to_t3_label
                }

                field_map = table_to_map_dict[table_num]

                extracted_data = {}
                for zip_file in zip_list:
                    with zipfile.ZipFile(zip_file) as myzip:
                        #remove .zip extension to get path within unzipped archive
                        zipfile_basename = os.path.basename(zip_file)
                        zipfile_nameroot = os.path.splitext(zipfile_basename)[0]
                        datafile_archive_path = zipfile_nameroot + '/fastqc_data.txt'
                        with myzip.open(datafile_archive_path) as myfile:
                            #the file is split into sections using this string as a
                            #delimiter; all required info is in the first section
                            raw_chunk = myfile.read().decode().split('>>END_MODULE')[0]

                    #split chunk on newlines, then split each line on tabs
                    parsed_data = [line.split('\t') for line in raw_chunk.splitlines()]
                    #ensure filename, which is used as a key for sub-dictionaries,
                    #is set once & only once; it's the first data field in the chunk,
                    #so the sub-dict is safely created before it's used to store data
                    key_unset = True
                    for (field_name, field_data) in parsed_data:
                        if key_unset and field_name == 'Filename':
                            key = field_data
                            extracted_data[key] = {}
                            key_unset = False
                        elif not key_unset and field_name in field_map:
                            extracted_data[key][ field_map[field_name] ] = field_data

                return extracted_data

            # the remaining methods produce single-level dictionaries, with keys
            # corresponding to table row names

            # the next two methods parse the outputs of Gue Su's custom metrics scripts

            def parse_aligned_metrics(metrics_file, table_num):
                fieldname_to_t2_label = {
                    "Total Read Count (R1 + R2)": "Total Read Count",
                    "Unique Read Pairs": "Unique Read Pairs (%)",
                    "Total Mapped Reads": "Total Mapped Reads (%)",
                    "Non-Mapped Reads": "Non-Mapped Reads (%)",
                    "Unique Mapped Reads": "Unique Mapped Reads (%)",
                    "Mapped Read Duplication": "Mapped Read Duplication (%)",
                    "Strand ratio (forward, reverse, reverse/forward of unique mapped)": "Strand Specificity (%)"
                }

                fieldname_to_t3_label = {
                    "QC-failed Read Count": "Filtered (read count)",
                    "Total Mapped Reads": "Total mapped reads (%)",
                    "Non-Mapped Reads": "Non-mapped reads (%)",
                    "Unique Mapped Reads": "Unique mapped reads (%)",
                    "Strand ratio (forward, reverse, reverse/forward of unique mapped)": "Strand specificity (%)"
                }

                table_to_map_dict = {
                    "table2": fieldname_to_t2_label,
                    "table3": fieldname_to_t3_label
                }

                field_map = table_to_map_dict[table_num]

                with open(metrics_file) as f:
                    parsed_data = [line.split('\t') for line in f.read().splitlines()]
                extracted_data = {}
                for parsed_line in parsed_data:
                    if len(parsed_line) < 2:
                        continue
                    field_name = parsed_line[0]
                    field_data = parsed_line[-1]
                    if field_name in field_map:
                        extracted_data[ field_map[field_name] ] = field_data

                return extracted_data

            def parse_unaligned_metrics(metrics_files):

                field_map = {
                    "Median Basecall Quality Score": "Median Phred Score",
                    "Bases with N": "Bases with N (%)",
                    "Bases with >= Q30": "Bases with >= Q30 (%)"
                }
                extracted_data = {}
                for metrics_file in metrics_files:
                    with open(metrics_file) as f:
                        parsed_data = [line.split('\t') for line in f.read().splitlines()]
                    file_data = {}
                    filenames = []
                    for parsed_line in parsed_data:
                        if 'cumulative' in parsed_line[0]:
                            filenames.append(os.path.basename(parsed_line[0].split()[-1]))

                        if len(parsed_line) < 2:
                            continue
                        field_name = parsed_line[0]
                        field_data = parsed_line[-1]
                        if field_name in field_map:
                            file_data[ field_map[field_name] ] = field_data
                    for fname in filenames:
                        extracted_data[fname] = file_data

                return extracted_data

            # the remaining methods parse various picard output files

            def parse_alignment_summary_metrics(alignment_summary_metrics):
                with open(alignment_summary_metrics) as f:
                    parsed_data = [line.split('\t') for line in f.read().splitlines()]
                extracted_data = {}
                for parsed_line in parsed_data:
                    if parsed_line[0] == 'PAIR':
                        extracted_data['PCT_READS_ALIGNED_IN_PAIRS'] = parsed_line[17]
                    elif parsed_line[0] == 'FIRST_OF_PAIR':
                        extracted_data['PF_MISMATCH_RATE_1'] = parsed_line[12]
                    elif parsed_line[0] == 'SECOND_OF_PAIR':
                        extracted_data['PF_MISMATCH_RATE_2'] = parsed_line[12]
                return extracted_data

            def indices_from_field_list(header_line, fields):
                field_to_index = {}
                for index, field in enumerate(header_line):
                    if field in fields:
                        field_to_index[field] = index
                        fields.remove(field)
                assert(len(fields) == 0)
                return field_to_index

            def parse_duplication_metrics(duplication_metrics):
                with open(duplication_metrics) as f:
                    raw_chunk = f.read().split('\n\n')[1]
                pct_dup = raw_chunk.splitlines()[2].split('\t')[8]
                return {'PERCENT_DUPLICATION': pct_dup}

            def parse_insert_size_metrics(insert_size_metrics):
                with open(insert_size_metrics) as f:
                    raw_chunk = f.read().split('\n\n')[1]
                insert_size = raw_chunk.splitlines()[2].split()[5]
                return {'MEAN_INSERT_SIZE': insert_size}

            def parse_hs_metrics(hs_metrics):
                with open(hs_metrics) as f:
                    raw_chunk = f.read().split('\n\n')[1]
                parsed_chunk = [line.split() for line in raw_chunk.splitlines()]
                header_line = parsed_chunk[1]
                data_line = parsed_chunk[2]
                fields = set(['PCT_USABLE_BASES_ON_TARGET', 'PCT_EXC_OFF_TARGET', 'MEAN_TARGET_COVERAGE', 'PCT_TARGET_BASES_20X'])
                field_to_index = indices_from_field_list(header_line, fields)
                extracted_data = {}
                for field in field_to_index:
                    extracted_data[field] = data_line[ field_to_index[field] ]
                return extracted_data
                
            def parse_rna_metrics(rna_metrics):
                with open(rna_metrics) as f:
                    raw_chunk = f.read().split('\n\n')[1]
                parsed_chunk = [line.split() for line in raw_chunk.splitlines()]
                header_line = parsed_chunk[1]
                data_line = parsed_chunk[2]
                fields = set([
                    'PCT_CODING_BASES', 'PCT_UTR_BASES', 'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES', 'PCT_CORRECT_STRAND_READS',
                    'MEDIAN_5PRIME_BIAS', 'MEDIAN_3PRIME_BIAS', 'MEDIAN_5PRIME_TO_3PRIME_BIAS'])
                field_to_index = indices_from_field_list(header_line, fields)
                extracted_data = {}
                for field in field_to_index:
                    extracted_data[field] = data_line[ field_to_index[field] ]
                return extracted_data

            def parse_flagstat(flagstat):
                with open(flagstat) as f:
                    return {'Filtered Read Count': f.readlines()[0].split()[2]}

            def parse_unaligned_rna_table(table_file):
                with open(table_file) as f:
                    parsed_data = [line.split(',') for line in f.read().splitlines()]
                for parsed_line in parsed_data:
                    if parsed_line[0] == 'Total Number of Reads':
                        total_reads = sum([int(num) for num in parsed_line[1:]])
                        return {'Total read (read count)': total_reads}

            def aggregate_dicts(md5_dict, fastqc_dict, *flat_dicts, unaligned_dict=None, rename_key=False):
                final_dict = md5_dict.copy()
                for key in final_dict:
                    alternate_key = key
                    if rename_key:
                        nameroot, nameext = os.path.splitext(key)
                        if nameext == '.cram':
                            alternate_key = nameroot + '.bam'
                    # merge the nested dicts, which are in the form {filename: {table_field: val, ...}}
                    final_dict[key].update(fastqc_dict[alternate_key])
                    if unaligned_dict:
                        final_dict[key].update(unaligned_dict[key])
                    # copy in any given flat dicts, in the form {table_field: val, ...}
                    # vals in flat dicts not keyed to filename are assumed to apply across all files
                    for flat_dict in flat_dicts:
                        final_dict[key].update(flat_dict)
                return final_dict

            def generate_table(table_rows, val_dict):
                table_list = []
                for row in table_rows:
                    row_list = [row]
                    for filename_key in sorted(val_dict):
                        row_list.append(val_dict[filename_key][row])
                    table_list.append(row_list)
                return table_list

            def generate_table1(file_args, string_arg_dict):
                # NOTE: 'Total Number of Reads' is used as a key in function parse_unaligned_rna_table
                # make sure to update that function if the name is changed below
                table_rows = ('Sample Name', 'File', 'MD5 File Checksum', 'Sequencing Platform',
                    'Sequencing Instrument', 'Sequencing Kit', 'Single or Paired End', 'Sequencing Type',
                    'Total Number of Reads', 'Median Phred Score', 'GC Content (%)', 'Bases with N (%)',
                    'Bases with >= Q30 (%)', 'Spike in Error Rate (%)', 'FastQC Plots Provided?')

                md5_dict = parse_md5sums(file_args.md5sum_file)
                fastqc_dict = parse_fastqc_zips(file_args.fastqc_zips, file_args.table)
                unaligned_dict = parse_unaligned_metrics(file_args.unaligned_metrics)

                table_dict = aggregate_dicts(md5_dict, fastqc_dict, string_arg_dict, unaligned_dict=unaligned_dict)

                table = generate_table(table_rows, table_dict)
                return table

            def generate_table2(file_args, string_arg_dict):
                table_rows = ('Sample Name', 'File', 'MD5 File Checksum', 'Sequencing Platform',
                    'Sequencing Instrument', 'Sequencing Kit', 'Single or Paired End', 'Source', 'Total DNA(ng)',
                    'Read Length', 'Total Read Count', 'Filtered Read Count', 'Unique Read Pairs (%)', 'Total Mapped Reads (%)',
                    'Non-Mapped Reads (%)', 'Unique Mapped Reads (%)', 'Mapped Read Duplication (%)', 'Strand Specificity (%)',
                    'Reference Genome', 'PCT_USABLE_BASES_ON_TARGET', 'PCT_EXC_OFF_TARGET', 'PERCENT_DUPLICATION',
                    'MEAN_TARGET_COVERAGE', 'PCT_TARGET_BASES_20X', 'PCT_READS_ALIGNED_IN_PAIRS', 'MEAN_INSERT_SIZE',
                    'PF_MISMATCH_RATE_1', 'PF_MISMATCH_RATE_2')

                md5_dict = parse_md5sums(file_args.md5sum_file)
                fastqc_dict = parse_fastqc_zips(file_args.fastqc_zips, file_args.table)
                aligned_dict = parse_aligned_metrics(file_args.aligned_metrics, file_args.table)
                alignment_summary_dict = parse_alignment_summary_metrics(args.alignment_summary_metrics)
                duplication_dict = parse_duplication_metrics(args.duplication_metrics)
                insert_size_dict = parse_insert_size_metrics(args.insert_size_metrics)
                hs_dict = parse_hs_metrics(args.hs_metrics)
                flagstat_dict = parse_flagstat(args.flagstat)

                table_dict = aggregate_dicts(md5_dict, fastqc_dict, string_arg_dict, aligned_dict,
                    alignment_summary_dict, duplication_dict, insert_size_dict, hs_dict, flagstat_dict, rename_key=True)

                table = generate_table(table_rows, table_dict)
                return table

            def generate_table3(file_args, string_arg_dict):
                table_rows = ('Sample Name', 'File', 'MD5 File Checksum', 'Sequencing Platform',
                    'Sequencing Instrument', 'Sequencing Kit', 'Single or Paired End', 'Source', 'Total RNA (ng)',
                    'RIN score', 'Read Length (nt)', 'Total read (read count)', 'Filtered (read count)', 'Total mapped reads (%)',
                    'Non-mapped reads (%)', 'Unique mapped reads (%)', 'Strand specificity (%)', 'Frequency Normalization Method',
                    'Annotation File', 'Reference Genome', 'PCT_CODING_BASES', 'PCT_UTR_BASES', 'PCT_INTRONIC_BASES', 'PCT_INTERGENIC_BASES',
                    'PCT_CORRECT_STRAND_READS', 'MEDIAN_5PRIME_BIAS', 'MEDIAN_3PRIME_BIAS', 'MEDIAN_5PRIME_TO_3PRIME_BIAS')

                md5_dict = parse_md5sums(file_args.md5sum_file)
                fastqc_dict = parse_fastqc_zips(file_args.fastqc_zips, file_args.table)
                aligned_dict = parse_aligned_metrics(file_args.aligned_metrics, file_args.table)
                rna_dict = parse_rna_metrics(args.rna_metrics)
                upstream_readcount_dict = parse_unaligned_rna_table(args.unaligned_rna_table)

                table_dict = aggregate_dicts(md5_dict, fastqc_dict, string_arg_dict, aligned_dict, rna_dict, upstream_readcount_dict)

                table = generate_table(table_rows, table_dict)
                return table

            if args.table == 'table1':
                table = generate_table1(args, args_to_table_fields)
            elif args.table == 'table2':
                table = generate_table2(args, args_to_table_fields)
            else:
                table = generate_table3(args, args_to_table_fields)

            with open(args.table_file_name, 'w+') as f:
                for line in table:
                    f.write(','.join([str(e) for e in line]) + '\n')


baseCommand: ['python', 'generate_tables.py']
inputs:
    table_file_name:
        type: string?
        inputBinding:
            prefix: "--table_file_name"
    table_num:
        type: string
        inputBinding:
            prefix: "--table"
    md5sum_file:
        type: File
        inputBinding:
            prefix: "--md5sum_file"
    fastqc_zips:
        type: File[]
        inputBinding:
            prefix: "--fastqc_zips"
    aligned_metrics:
        type: File?
        inputBinding:
            prefix: "--aligned_metrics"
    unaligned_metrics:
        type: File[]?
        inputBinding:
            prefix: "--unaligned_metrics"
    alignment_summary_metrics:
        type: File?
        inputBinding:
            prefix: "--alignment_summary_metrics"
    duplication_metrics:
        type: File?
        inputBinding:
            prefix: "--duplication_metrics"
    insert_size_metrics:
        type: File?
        inputBinding:
            prefix: "--insert_size_metrics"
    hs_metrics:
        type: File?
        inputBinding:
            prefix: "--hs_metrics"
    rna_metrics:
        type: File?
        inputBinding:
            prefix: "--rna_metrics"
    flagstat:
        type: File?
        inputBinding:
            prefix: "--flagstat"
    unaligned_rna_table:
        type: File?
        inputBinding:
            prefix: "--unaligned_rna_table"

    sequencing_platform:
        type: string?
        inputBinding:
            prefix: "--sequencing_platform"
    sequencing_instrument:
        type: string?
        inputBinding:
            prefix: "--sequencing_instrument"
    sequencing_kit:
        type: string?
        inputBinding:
            prefix: "--sequencing_kit"
    sequencing_type:
        type: string?
        inputBinding:
            prefix: "--sequencing_type"
    single_or_paired_end:
        type: string?
        inputBinding:
            prefix: "--single_or_paired_end"
    spike_in_error_rate:
        type: string?
        inputBinding:
            prefix: "--spike_in_error_rate"
    source:
        type: string?
        inputBinding:
            prefix: "--source"
    total_DNA:
        type: string?
        inputBinding:
            prefix: "--total_DNA"
    reference_genome:
        type: string?
        inputBinding:
            prefix: "--reference_genome"
    total_RNA:
        type: string?
        inputBinding:
            prefix: "--total_RNA"
    RIN_score:
        type: string?
        inputBinding:
            prefix: "--RIN_score"
    freq_normalization_method:
        type: string?
        inputBinding:
            prefix: "--freq_normalization_method"
    annotation_file:
        type: string?
        inputBinding:
            prefix: "--annotation_file"
    sample_name:
        type: string?
        inputBinding:
            prefix: "--sample_name"
outputs:
    table:
        type: File
        outputBinding:
            glob: $(inputs.table_file_name)
