cwlVersion: v1.2
class: CommandLineTool

baseCommand: python

inputs:
  script:
    type: File
    default: 
      class: File
      location: "../../modules/ranking/select_and_run.py"
    inputBinding:
      position: -1
  queries:
    type: File
    inputBinding:
      position: 1
  expression_matrix:
    type: File
    inputBinding:
      position: 2
  metadata_matrix:
    type: File
    inputBinding:
      position: 3
  method:
    type: string
    inputBinding:
      prefix: --method

outputs:
  rankings:
    type:
      type: array
      items: [File]
    outputBinding:
      glob: "*_deseq.csv"

arguments:
  - position: 4
    valueFrom: $(runtime.outdir)
