cwlVersion: v1.0
class: CommandLineTool
label: "add population allele frequencies to a vcf"

arguments: ["/opt/hall-lab/python-2.7.15/bin/svtools", "afreq", "$(inputs.vcf)", { shellQuote: false, valueFrom: "|" },
    "/opt/hall-lab/python-2.7.15/bin/svtools", "vcftobedpe", "-i", "stdin", { shellQuote: false, valueFrom: "|" },
    "/opt/hall-lab/python-2.7.15/bin/svtools", "varlookup", "-d", "200", "-c", "POPFREQ", "-a", "stdin", "-b", "$(inputs.sv_db)", { shellQuote: false, valueFrom: "|" },
    "/opt/hall-lab/python-2.7.15/bin/svtools", "bedpetovcf"
    ]
stdout: "$(inputs.cohort_name)"
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "halllab/svtools:v0.4.0-2dc7f91" 
    - class: ResourceRequirement
      ramMin: 5000
      coresMin: 1

inputs:
    vcf:
        type: File
        doc: "vcf to add populationa allele frequence to"
    sv_db:
        type: File
        doc: "bed file containing allele frequencies for a population"
    cohort_name:
        type: string?
        default: "merged_annotated.vcf"

outputs:
    merged_annotated_vcf:
        type: stdout
