#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: ExpressionTool
label: "Renamer"
doc: "Logically renames a file, but not all workflow engines currently support it.  See also staged_rename.cwl"
requirements:
    - class: InlineJavascriptRequirement
inputs:
    original:
        type: File
    name:
        type: string
outputs:
    replacement:
        type: File
expression: |
    ${
        var f = inputs.original;
        f.basename = inputs.name;
        return {'replacement': f};
    }
