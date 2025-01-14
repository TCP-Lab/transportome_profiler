cwlVersion: v1.2
class: CommandLineTool

baseCommand: generanker

inputs:
  id_column:
    type: string
    inputBinding:
      prefix: --id-col
  case_matrix:
    type: File
    inputBinding:
      position: 1
  control_matrix:
    type: File
    inputBinding:
      position: 2
  metric:
    type: string
    inputBinding:
      position: 10

arguments:
  - prefix: "--output-file"
    valueFrom: "output.csv"

outputs:
  rank_matrix:
    type: File
    outputBinding:
      glob: "output.csv"
