cwlVersion: v1.0
class: CommandLineTool
label: "Concordance tool"
baseCommand: ["python", "/opt/concordance/newConcordance.py"]
#commented out bc currently cannot do
#requirements:
#    - class: DockerRequirement
#      dockerPull: "mneveau/docker-cle:concordance_addition"
inputs:
    bam_1:
        type: File
        inputBinding:
            position: 2
        secondaryFiles: [.bai]
    bam_2:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [.bai]
    reference:
        type: File
        inputBinding:
            position: 4
        secondaryFiles: [.fai, "^.dict"]
    snp:
        type: File
        inputBinding:
            position: 5
    output_file_name:
        type: string?
        inputBinding:
            position: 6
	default: "concordance_report.txt" 
outputs:
    output_file:
        type: File
        outputBinding:
            glob: "*.txt"