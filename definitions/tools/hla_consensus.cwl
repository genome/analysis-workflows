#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: "Script to create consensus from optitype and clinical HLA typing"
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "python:3.7.4-slim-buster"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'hla_consensus.py'
        entry: |
                        #This script produces 2-4 files depending on inputs and their contents
                        #All are packaged together into a folder called hla_calls for convenience
                        #optitype_calls.txt is always produced, and is essentially a copy of optitype's output
                        #consensus_calls.txt is also always produced; if no clinical calls are provided, this
                        #file is identical to optitype_calls.txt. If clinical calls are provided, they are
                        #reproduced in clinical_calls.txt. If the clinical calls exactly match the optitype calls*, 
                        #all 3 files described so far will contain the same information, but are not guaranteed to
                        #be exactly the same (text ordering may differ, depending on the order calls are given in the input). 
                        #If the clinical calls and optitype calls do not match, mismatched_calls.txt is then produced;
                        #each line represents a gene. See below (section 'write out call files') for more mismatch details.
                        #NOTE: optitype only produces MHC class I calls, 

                        #optitype input format:
                        #HLA-[A,B,C, ...]*xx:xx
                        #clinical input format:
                        # HLA-C*12:02/HLA-C*12:228/HLA-C*12:243 (everything first slash and after is optional)

                        import sys, os
                        from collections import defaultdict

                        ####################################
                        ### helper methods for later use ###
                        ####################################

                        #helper method that takes in the decomposed version of an hla
                        #string and returns the full delimited string
                        def build_hla_str(gene, allele_group, spec_allele):
                            return gene + "*" + allele_group + ":" + spec_allele

                        #helper method that takes in a full hla string, like HLA-X*01:02:03:04,
                        #and splits it into the gene name (HLA-X), allele group (01), and the
                        #specific allele (02), dropping any fields beyond this, because downstream
                        #tool do not support these fields
                        def split_hla_str(full_hla_str):
                            gene_name, raw_allele_fields = full_hla_str.split('*')
                            split_allele_fields = raw_allele_fields.split(":")
                            allele_group_name = split_allele_fields[0]
                            specific_allele_name = split_allele_fields[1]
                            return (gene_name, allele_group_name, specific_allele_name)

                        #helper method that creates a mismatch file only if any have been found in the tree,
                        #and inserts a header upon initially creating the file. Params:
                        #previously_written- true if the file has already been created; used control header creation
                        #mismatches- dictionary with sources as keys and a list of alleles called only by that
                        #            source as values
                        #returns true if the file was or has ever been written to, false otherwise
                        def write_mismatch(previously_written, mismatches):
                            if (not(mismatches['optitype'] or mismatches['clinical'])):
                                #In this case, both arrays are empty, so there's no mismatch to write
                                #function has not changed the file state, so return the unmodified flag
                                return previously_written

                            with open("hla_calls/mismatched_calls.txt", "a") as m_c:
                                if not previously_written:
                                    #add header if this is the first time writing to the file
                                    m_c.write("optitype_calls\tclinical_calls\n")
                                #write the mismatches to the file
                                m_c.write( ",".join(mismatches['optitype']) + "\t" + ",".join(mismatches['clinical']) + "\n" )

                            return True

                        ########################################
                        ### parse args from the command line ###
                        ########################################

                        clinical_exists = len(sys.argv) > 2 

                        optitype_calls = sys.argv[1].split(",")

                        if clinical_exists:
                            raw_clinical_calls = sys.argv[2].split(",")
                            #Each clinical call may be a single high confidence call,
                            #or a list of uncertain calls separated by slashes
                            hc_clinical_calls = []
                            u_clinical_calls = []
                            for call in raw_clinical_calls:
                                if "/" in call:
                                    u_clinical_calls.append(call)
                                else:
                                    hc_clinical_calls.append(call)

                        ################################################################
                        ### Load HLA types into data structure for consensus calling ###
                        ################################################################

                        #Create a basic tree out of dictionaries to hold the data from all callers;
                        #top level keys will be genes, pointing to a nested dictionary with
                        #allele groups as keys, pointing to a final dictionary with specific alleles
                        #as keys and a set containing call sources as values
                        # ex: optitype calls HLA-A*01:02 -> {HLA-A: {01: {02: {optitype}}}}

                        #defaultdict constructor requires a callable; however, it returns an object
                        #lambdas create a callable that allows for nested defaultdicts
                        hla_calls = defaultdict( lambda: defaultdict( lambda: defaultdict(set) ) )

                        for call in optitype_calls:
                            gene, allele_group, spec_allele = split_hla_str(call)

                            #records this call in the tree and tags it as coming from optitype
                            #the tag is added to a set, so any duplicates (such as from an individual homozygous
                            #for a given gene) are collapsed into a single entry
                            hla_calls[gene][allele_group][spec_allele].add('optitype')

                        if clinical_exists:
                            for call in hc_clinical_calls:
                                gene, allele_group, spec_allele = split_hla_str(call)

                                #Case 1: this $call was also called by optitype, so add to the 
                                #record indicating that clinical data supports this call
                                #Case 2: this call is unique to the clinical data; create a record 
                                #and indicate that only clinical data supports this call
                                hla_calls[gene][allele_group][spec_allele].add('clinical')

                            for multi_call in u_clinical_calls:
                                calls = multi_call.split("/")
                                multi_consensus = set()
                                for call in calls:
                                    gene, allele_group, spec_allele = split_hla_str(call)

                                    #check if this call already exists in the tree, which will be treated as
                                    #evidence that this call is the correct call out of the current group of
                                    #uncertain calls ($multi_call)

                                    #TODO this may be biased towards creating a homozygous consensus:
                                    #since high confidence clinical calls are evaluated before this, one of 
                                    #these calls could be used as evidence when resolving the uncertain call
                                    #using the current method. Is this desirable? Should we only use optitype calls
                                    #when resolving uncertain clinical calls?
                                    #Example: clinical calls 01:02 and 01:02/01:03/01:04
                                    if hla_calls[gene][allele_group][spec_allele]:
                                        #add as a tuple to avoid re-splitting later
                                        multi_consensus.add( (gene, allele_group, spec_allele) )

                                #if one and only one of the calls from the uncertain group was already in the tree,
                                #that is treated as evidence that this particular call was the correct one. It will
                                #be accepted and entered into the tree, while the other calls will be discarded
                                if len(multi_consensus) == 1:
                                    accpt_call = multi_consensus.pop()
                                    hla_calls[accpt_call[0]][accpt_call[1]][accpt_call[2]].add('clinical')
                                #otherwise, all uncertain calls from the group will be added to the tree; this means
                                #they will be added to the consensus superset (since their validity cannot be disproven),
                                #and also used to construct the mismatch file
                                else:
                                    for call in calls:
                                        gene, allele_group, spec_allele = split_hla_str(call)
                                        hla_calls[gene][allele_group][spec_allele].add('clinical')

                        ##############################################
                        ### write out caller files for convenience ###
                        ##############################################

                        os.mkdir("hla_calls")

                        #Create an exact copy of optitype calls, to be bundled with other relevant
                        #files for convenience/later review. Always generated,
                        with open("hla_calls/optitype_calls.txt", "w") as o_c:
                            o_c.write( ",".join(optitype_calls) )

                        #Create an exact copy of clinical calls, if they exist, to be bundled with 
                        #other relevant files for convenience/later review.
                        if clinical_exists:
                            with open("hla_calls/clinical_calls.txt", "w") as c_c:
                                c_c.write( ",".join(raw_clinical_calls) )

                        #########################################################
                        ### Generate consensus (superset if callers disagree) ###
                        #########################################################

                        #A consensus file is always generated to be passed on to pvacseq. If there are
                        #no clinical calls, this file is the same as optitype_calls.txt. If there are, walk
                        #through the tree and emit everything present as the consensus. If there is a true
                        #consensus, each class I gene (corresponding to the top level keys of the tree) will have
                        #at most 2 leaves (1 in the case of a homozygote, or in the rare case that both optitype
                        #and clinical data only called one allele for this gene), where each leaf represents
                        #a specific allele call supported by both sources. If there is no true consensus, there
                        #may be more than 2 leaves per class I gene, and individual leaves may only be supported by
                        #1 of the 2 sources. These leaves will still be added to the consensus to form a superset,
                        #since there is not enough evidence to discard them, but they will also be added to a
                        #mismatch file, which presents side by side lists of the differing alleles called by each
                        #source, with one gene per line. Note that optitype only makes class I predictions, so any
                        #class II predictions from the clinical data are always added to the consensus and never
                        #to the mismatch file
                        if not clinical_exists:
                            with open("hla_calls/consensus_calls.txt", "w") as c_c:
                                c_c.write( ",".join(optitype_calls) )
                        else:
                            consensus_calls = []
                            mismatch_written = False
                            for gene in hla_calls:
                                mismatches = {'optitype': [], 'clinical': []}
                                for allele_group in hla_calls[gene]:
                                    for spec_allele in hla_calls[gene][allele_group]:
                                        callers = hla_calls[gene][allele_group][spec_allele]

                                        #if any uncertain calls were resolved to a single call based on prior
                                        #evidence, the discarded calls will have been visited but not tagged,
                                        #resulting in leaves with empty sets; these can be ignored
                                        if callers:
                                            #there are now only 3 possibilities for the contents of $callers:
                                            #{optitype, clinical}, {optitype}, {clinical}
                                            #all will be added to the consensus, possibly creating a superset
                                            #those with only 1 caller represent mismatches between the 2
                                            consensus_calls.append( build_hla_str(gene, allele_group, spec_allele) )
                                            if len(callers) == 1:
                                                mismatches[callers.pop()].append( build_hla_str(gene, allele_group, spec_allele) )

                                mismatch_written = write_mismatch(mismatch_written, mismatches)

                            with open("hla_calls/consensus_calls.txt", "w") as c_c:
                                c_c.write( ",".join(consensus_calls) )

baseCommand: ['python', 'hla_consensus.py']
inputs:
    optitype_calls:
        type: string[]
        inputBinding:
            position: 1
            itemSeparator: ','
            separate: false
    clinical_calls:
        type: string[]?
        inputBinding:
            position: 2
            itemSeparator: ','
            separate: false
outputs:
    consensus_alleles:
        type: string[]
        outputBinding:
            glob: hla_calls/consensus_calls.txt
            loadContents: true
            outputEval: $(self[0].contents.split(","))
    hla_call_files:
        type: Directory
        outputBinding:
            glob: hla_calls
