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
                        import sys, os
    
                        clinical_exists = len(sys.argv) > 2 
                        consensus = []
                        mismatched = False
    
                        optitype_calls = sys.argv[1].split(",")
                        if clinical_exists:
                            clinical_calls = sys.argv[2].split(",")
    
                            a_alleles = [ set(), set() ]
                            b_alleles = [ set(), set() ]
                            c_alleles = [ set(), set() ]
    
                        os.mkdir("hla_calls")
    
                        with open("hla_calls/optitype_calls.txt", "w") as o_c:
                            o_c.write( ",".join(optitype_calls) )
    
                            for call in optitype_calls:
    
                                if clinical_exists:
                                    if "HLA-A" in call:
                                        a_alleles[0].add(call)
                                    elif "HLA-B" in call:
                                        b_alleles[0].add(call)
                                    elif "HLA-C" in call:
                                        c_alleles[0].add(call)
    
                        if clinical_exists:
                            with open("hla_calls/clinical_calls.txt", "w") as c_c:
                                c_c.write( ",".join(clinical_calls) )
    
                                for call in clinical_calls:
    
                                    if "HLA-A" in call:
                                        if call in a_alleles[0]:
                                            a_alleles[0].remove(call)
                                        else:
                                            mismatched = True
                                            a_alleles[1].add(call)
                                    elif "HLA-B" in call:
                                        if call in b_alleles[0]:
                                            b_alleles[0].remove(call)
                                        else:
                                            mismatched = True
                                            b_alleles[1].add(call)
                                    elif "HLA-C" in call:
                                        if call in c_alleles[0]:
                                            c_alleles[0].remove(call)
                                        else:
                                            mismatched = True
                                            c_alleles[1].add(call)
    
                                    consensus.append(call)
                        else:
                            consensus = optitype_calls
    
                        with open("hla_calls/consensus_calls.txt", "w") as c_c:
                            c_c.write( ",".join(consensus) )
    
                        if mismatched:
                            with open("hla_calls/mismatched_calls.txt", "w") as m_c:
                                m_c.write( "optitype_calls\tclinical_calls\n" )
                                m_c.write( ",".join(a_alleles[0]) + "\t" + ",".join(a_alleles[1]) + "\n" )
                                m_c.write( ",".join(b_alleles[0]) + "\t" + ",".join(b_alleles[1]) + "\n" )
                                m_c.write( ",".join(c_alleles[0]) + "\t" + ",".join(c_alleles[1]) )            


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
