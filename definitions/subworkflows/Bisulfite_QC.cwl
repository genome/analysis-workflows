#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Bisulfite QC tools"
requirements:
    - class: SubworkflowFeatureRequirement

inputs:
      vcf:
          type: File
      bam:
          type: File
      reference:
          type: string
      QCannotation:
           type: File
      output_dir:
           type: string?
           default: bisulfite_QC

outputs:
    QC_directory:
        type: Directory
        outputSource: gather_to_sub_directory/gathered_directory

steps:
      Tool_bisulfiteconversion:
          run: ../tools/Tool_bisulfiteconversion.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
              [baseconvertion, readconvertion, CpHretention, CpGretention]

      Tool_mappingsummary:
          run: ../tools/Tool_mappingsummary.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
             [Strandtable, Mappingquality]

      Tool_CpGretentiondistribution:
          run: ../tools/Tool_CpGretentiondistribution.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
             [CpGRetentionDist]

      Tool_coveragestats:
          run: ../tools/Tool_coveragestats.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
             [BgaBed, covdist, BgaBeddup, Dupreport, CpGbed, covdistcpg, cpgdist]

      gather_to_sub_directory:
          run: ../tools/gather_to_sub_directory.cwl
          in:
            outdir: output_dir
            files: [Tool_bisulfiteconversion/baseconvertion, Tool_bisulfiteconversion/readconvertion,
                  Tool_bisulfiteconversion/CpHretention, Tool_bisulfiteconversion/CpGretention,
                  Tool_mappingsummary/Strandtable, Tool_mappingsummary/Mappingquality,
                  Tool_CpGretentiondistribution/CpGRetentionDist,
                  Tool_coveragestats/BgaBed, Tool_coveragestats/covdist, Tool_coveragestats/BgaBeddup, Tool_coveragestats/Dupreport,
                  Tool_coveragestats/CpGbed, Tool_coveragestats/covdistcpg, Tool_coveragestats/cpgdist]
          out:
             [gathered_directory]
