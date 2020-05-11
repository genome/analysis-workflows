class: ExpressionTool
# takes a file and an array of files and prepends the file to the array
cwlVersion: v1.0
inputs:
  file:
    type: File
  array:
    type: File[]
outputs:
    file_array: File[]
requirements:
    InlineJavascriptRequirement: {}
expression: |
   ${
      var afile = inputs.file;
      var file_array = inputs.array;
      file_array.unshift(afile);
      return {"file_array": file_array}
    }
