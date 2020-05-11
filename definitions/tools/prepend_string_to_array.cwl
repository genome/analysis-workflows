class: ExpressionTool
# takes a string and an array of strings and prepends the string to the array
cwlVersion: v1.0
inputs:
  string:
    type: string
  array:
    type: string[]
outputs:
    string_array: string[]
requirements:
    InlineJavascriptRequirement: {}
expression: |
   ${
      var astring = inputs.string;
      var string_array = inputs.array;
      string_array.unshift(astring);
      return {"string_array": string_array}
    }
