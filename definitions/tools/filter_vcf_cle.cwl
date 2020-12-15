#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Filter variants from the VCF annotated with CLE-validated somatic VCF"
baseCommand: ["/usr/bin/perl", "cle_filter.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/cle:v1.4.2"
    - class: ResourceRequirement
      ramMin: 4000
    - class: StepInputExpressionRequirement
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'cle_filter.pl'
        entry: |
            #!/usr/bin/perl
            
            # This script intends to prepare a somatic VCF input file 
            # for pVACseq analysis with clinical trials.
            # For example, it updates the FILTER status of the variant whose FILTER status is PASS
            # but, if that variant is not validated by the CLE somatic pipeline.
            # In sum, the goal of this script is
            # to feed only CLE-validated somatic variants to pVACseq analysis.
            # The final VCF output is formatted as pVACseq requires.
            
            # It is free software; you can redistribute it
            # and/or modify it under the same terms as Perl itself.
            # This software is provided "as is" without warranty of any kind.
            # See the The MIT License (MIT) for more details.
            
            
            use strict;
            use warnings;
            
            use feature qw(say);
            
            # parses user's arguments: input_vcf_path, filter_flag
            my ($input_vcf_path, $filter_flag) = @ARGV;
            
            # defines the default for "input_vcf_path", "output_dir", "filter_flag" arguments
            # explicitly declared here if any
            die("Missing argument: input_vcf_path required") unless defined $input_vcf_path;
            $filter_flag = 0 unless defined $filter_flag;
            
            
            my $tag = 'PREVIOUSLY_DISCOVERED';
            my $filter = 'PASS_NONCLE';
            
            # creates file handlers
            my $fh;
            if($input_vcf_path =~ /.gz$/) 
            {
                open($fh, "gunzip -c $input_vcf_path |")
                    or die("couldn't open $input_vcf_path to read");
            }
            else
            {
                open($fh, $input_vcf_path)
                    or die("couldn't open $input_vcf_path to read");
            }
            
            open (my $fout, ">", "annotated_filtered.vcf")
                or die "couldn't open annotated_filtered.vcf to write";
            
            
            while (<$fh>)
            {
                chomp;
                
                if (/^#/)
                {
                    say $fout $_;
                }
                elsif (/^#CHROM/)
                {
                    if ($filter_flag)
                    {
                        say $fout '##FILTER=<ID=' . $filter . ', Description="any variant of PASS filter status, but it has not been validated by CLE somatic pipeline">';
                    }
                    
                    # writes the header
                    say $fout $_;
                }
                else
                {
                    my @columns = split /\t/, $_;
                    
                    # the following source code was written by the current rules for CLE clinical cases
                    if ($filter_flag)
                    {
                        # from the INFO column
                        if ($columns[7] =~ /\;$tag/)
                        {
                            # writes a variant validated by the CLE somatic pipeline without any FILTER status change
                            say $fout $_;
                        }
                        else
                        {
                            # updates the filter status for any PASS variant calls that haven't been validated by CLE somatic pipeline
                            if ($columns[6] eq 'PASS')
                            {
                                $columns[6] = $filter;
                            }
                            elsif ($columns[6] =~ /PASS/ && $columns[6] =~ /\;/)
                            {
                                # in case of semicolon-separated multiple filter status
                                my @filters = split /\;/, $columns[6];
                                
                                for (my $i=0; $i<@filters; $i++)
                                {
                                    if ($filters[$i] eq 'PASS')
                                    {
                                        $filters[$i] = $filter;
                                    }
                                }
                                
                                # updates the FILTER column
                                $columns[6] = join ';', @filters;
                            }
                            #else
                            #{
                            #    # nothing done for nonPASS variants
                            #}
                            
                            say $fout join("\t", @columns);
                        }
                    }
                    else
                    {
                        # does nothing for variant filtration
                        say $fout $_;
                    }
                }
            }
            
            close($fh);
            close($fout);

inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    filter:
        type: boolean
        inputBinding:
            position: 2
            valueFrom: |
                ${
                  if(inputs.filter){
                    return "1";
                  } else {
                    return "0";
                  }
                }
outputs:
    cle_filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated_filtered.vcf"
