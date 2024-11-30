#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "VariantFiltration (GATK 4.1.8.1)"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx4g", "VariantFiltration"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
arguments:
    ["-O", { valueFrom: $(runtime.outdir)/$(inputs.output_vcf_basename).vcf.gz }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    vcf:
        type: File
        inputBinding:
            prefix: "--variant"
            position: 2
        secondaryFiles: [.tbi]
    filters:
        type: string[]
        inputBinding:
            position: 3
            valueFrom: |
              ${
                var results = []
                for(var i=0; i<self.length; i++){
                  var [filter, name] = self[i].split(";")
                  results.push("-filter")
                  results.push(filter)
                  results.push("--filter-name")
                  results.push(name)
                }
                return results
              }
        doc: "input array of strings with filter expression and filter name, split by a ';', Examples: 'QD<2.0;QD2', 'QUAL<30.0;QUAL30', 'SOR>3.0;SOR3'"
    output_vcf_basename:
        type: string?
        default: select_variants
outputs:
    filtered_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputBinding:
            glob: $(inputs.output_vcf_basename).vcf.gz
