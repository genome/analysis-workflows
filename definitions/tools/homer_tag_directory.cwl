#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Creating the homer tag directory'
doc: |
  homer annotation data is too large to include in the docker image, so it 
  requires a directory to be mounted at /opt/homerdata/ containing
    - config.txt - homer configuration file with directories pointing to paths like "data/accession"
    - data - folder containing homer annotation data files
  at WUSTL, this can be provided by providing the following in an analysis-project configuration:
    docker_volumes: "/gscmnt/gc2560/core/annotation_data/homer:/opt/homerdata"
  or outside the pipelines:
    LSF_DOCKER_VOLUMES="$LSF_DOCKER_VOLUMES /gscmnt/gc2560/core/annotation_data/homer:/opt/homerdata"
  for compute1, instead use the path `/storage1/fs1/bga/Active/gmsroot/gc2560/core/annotation_data/homer`
baseCommand: ["/bin/bash", "homer_tag_directory.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/homer:4.11.1"
    - class: ResourceRequirement
      ramMin: 32000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'homer_tag_directory.sh'
        entry: |
                set -o pipefail
                set -o errexit

                name="homer_tag_directory"
                bam=$1
                mkdir -p $name
                echo "creating tagDir"
                /opt/homer/bin/makeTagDirectory $name $bam

inputs:
    sam:
        type: File
        inputBinding:
            position: 1
outputs:
    tag_directory:
        type: Directory
        outputBinding:
            glob: "homer_tag_directory"
