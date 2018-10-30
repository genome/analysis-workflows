#! /usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cp']

requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"

arguments: ["-t"]

inputs:
    output_dir:
        type: string
        inputBinding:
            position: 1
    all_files:
        type: File[] 
        inputBinding:
            position: 2
outputs:
    gathered_files:
        type:
            type: array
            items: string
        outputBinding:
            outputEval: ${
                            var file_paths = [];
                            var new_path = inputs.output_dir;
                            if (new_path.slice(-1) != "/") {
                                new_path += "/";
                            }
                            var arrLen = inputs.all_files.length;
                            for(var i = 0; i < arrLen; i++) {
                                file_paths.push(new_path + inputs.all_files[i].basename);
                            }
                            return file_paths;
                        }
