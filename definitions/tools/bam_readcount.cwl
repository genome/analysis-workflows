#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bam-readcount"

baseCommand: ["/usr/bin/python", "bam_readcount_helper.py"]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/bam_readcount_helper-cwl:1.1.1"
    - class: ResourceRequirement
      ramMin: 16000
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'bam_readcount_helper.py'
        entry: |
            #!/usr/bin/env python

            import sys
            import os
            from cyvcf2 import VCF
            import tempfile
            import csv
            from subprocess import Popen, PIPE

            def generate_region_list(hash):
                fh = tempfile.NamedTemporaryFile('w', delete=False)
                writer = csv.writer(fh, delimiter='\t')
                for chr, positions in hash.items():
                    for pos in sorted(positions.keys()):
                        writer.writerow([chr, pos, pos])
                fh.close()
                return fh.name

            def filter_sites_in_hash(region_list, bam_file, ref_fasta, prefixed_sample, output_dir, insertion_centric, map_qual, base_qual):
                bam_readcount_cmd = ['/usr/bin/bam-readcount', '-f', ref_fasta, '-l', region_list, '-w', '0', '-b', str(base_qual), '-q', str(map_qual)]
                if insertion_centric:
                    bam_readcount_cmd.append('-i')
                    output_file = os.path.join(output_dir, prefixed_sample + '_bam_readcount_indel.tsv')
                else:
                    output_file = os.path.join(output_dir, prefixed_sample + '_bam_readcount_snv.tsv')
                bam_readcount_cmd.append(bam_file)
                execution = Popen(bam_readcount_cmd, stdout=PIPE, stderr=PIPE)
                stdout, stderr = execution.communicate()
                if execution.returncode == 0:
                    with open(output_file, 'wb') as output_fh:
                        output_fh.write(stdout)
                else:
                    sys.exit(stderr)

            #initializing these with default values
            min_base_qual = 20
            min_mapping_qual = 0

            if len(sys.argv) == 7:
                (script_name, vcf_filename, sample, ref_fasta, bam_file, prefix, output_dir)= sys.argv
            elif len(sys.argv) == 8:
                (script_name, vcf_filename, sample, ref_fasta, bam_file, prefix, output_dir, min_base_qual)= sys.argv
            elif len(sys.argv) == 9: #elif instead of else for explicit safety
                (script_name, vcf_filename, sample, ref_fasta, bam_file, prefix, output_dir, min_base_qual, min_mapping_qual)= sys.argv

            if prefix == 'NOPREFIX':
                prefixed_sample = sample
            else:
                prefixed_sample = '_'.join([prefix, sample])

            vcf_file = VCF(vcf_filename)
            sample_index = vcf_file.samples.index(sample)

            rc_for_indel = {}
            rc_for_snp   = {}
            for variant in vcf_file:
                ref = variant.REF
                chr = variant.CHROM
                start = variant.start
                end = variant.end
                pos = variant.POS
                for var in  variant.ALT:
                    if len(ref) > 1 or len(var) > 1:
                        #it's an indel or mnp
                        if len(ref) == len(var) or (len(ref) > 1 and len(var) > 1):
                            sys.stderr.write("Complex variant or MNP will be skipped: %s\t%s\t%s\t%s\n" % (chr, pos, ref , var))
                            continue
                        elif len(ref) > len(var):
                            #it's a deletion
                            pos += 1
                            unmodified_ref = ref
                            ref = unmodified_ref[1]
                            var = "-%s" % unmodified_ref[1:]
                        else:
                            #it's an insertion
                            var = "+%s" % var[1:]
                        if chr not in rc_for_indel:
                            rc_for_indel[chr] = {}
                        if pos not in rc_for_indel[chr]:
                            rc_for_indel[chr][pos] = {}
                        if ref not in rc_for_indel[chr][pos]:
                            rc_for_indel[chr][pos][ref] = {}
                        rc_for_indel[chr][pos][ref] = variant
                    else:
                        #it's a SNP
                        if chr not in rc_for_snp:
                            rc_for_snp[chr] = {}
                        if pos not in rc_for_snp[chr]:
                            rc_for_snp[chr][pos] = {}
                        if ref not in rc_for_snp[chr][pos]:
                            rc_for_snp[chr][pos][ref] = {}
                        rc_for_snp[chr][pos][ref] = variant

            if len(rc_for_snp.keys()) > 0:
                region_file = generate_region_list(rc_for_snp)
                filter_sites_in_hash(region_file, bam_file, ref_fasta, prefixed_sample, output_dir, False, min_mapping_qual, min_base_qual)
            else:
                output_file = os.path.join(output_dir, prefixed_sample + '_bam_readcount_snv.tsv')
                open(output_file, 'w').close()

            if len(rc_for_indel.keys()) > 0:
                region_file = generate_region_list(rc_for_indel)
                filter_sites_in_hash(region_file, bam_file, ref_fasta, prefixed_sample, output_dir, True, min_mapping_qual, min_base_qual)
            else:
                output_file = os.path.join(output_dir, prefixed_sample + '_bam_readcount_indel.tsv')
                open(output_file, 'w').close()

arguments: [
    { valueFrom: $(runtime.outdir), position: -3 }
]
stdout: $(inputs.sample)_bam_readcount.tsv
inputs:
    vcf:
        type: File
        inputBinding:
            position: -8
    sample:
        type: string
        inputBinding:
            position: -7
    reference_fasta:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            position: -6
    bam:
        type: File
        inputBinding:
            position: -5
        secondaryFiles: [.bai]
    prefix:
        type: string?
        default: 'NOPREFIX'
        inputBinding:
            position: -4
    min_base_quality:
        type: int?
        default: 20
        inputBinding:
            position: -2
    min_mapping_quality:
        type: int?
        default: 0
        inputBinding:
            position: -1
outputs:
    snv_bam_readcount_tsv:
        type: File
        outputBinding:
            glob: |
                    ${
                        var name = "_bam_readcount_snv.tsv";
                        if (inputs.prefix.equals("NOPREFIX")) {
                            name = inputs.sample + name;
                        }
                        else {
                            name = inputs.prefix + "_" + inputs.sample + name;
                        }
                        return name;
                    }
    indel_bam_readcount_tsv:
        type: File
        outputBinding:
            glob: |
                    ${
                        var name = "_bam_readcount_indel.tsv";
                        if (inputs.prefix.equals("NOPREFIX")) {
                            name = inputs.sample + name;
                        }
                        else {
                            name = inputs.prefix + "_" + inputs.sample + name;
                        }
                        return name;
                    }
            
