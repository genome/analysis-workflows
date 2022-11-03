#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: "Rename file paths with identical basenames so that output staging works"
requirements:
    - class: DockerRequirement
      dockerPull: "python:3.7.4-slim-buster"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'stage_for_renaming.py'
        entry: |
            import sys, os, shutil

            file_list = sys.argv[1:]

            for num, f in enumerate(file_list):
                original_name = os.path.basename(f)
                name_root, name_ext = os.path.splitext(original_name)
                new_name = name_root + '.' + str(num+1) + name_ext
                shutil.copy(f, new_name)

baseCommand: ['python', 'stage_for_renaming.py']
inputs:
    scatter_files:
        type: File[]
        inputBinding:
            position: 1
outputs:
    unique_files:
        type: File[]
        outputBinding:
            glob: "$(inputs.scatter_files[0].nameroot)*"
