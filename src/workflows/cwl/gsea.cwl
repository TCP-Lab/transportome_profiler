cwlVersion: v1.2
class: CommandLineTool

basecommand: python

inputs:
  script:
    type: File
    default:
      class: File
      location: ../../modules/run_gsea.R
    inputBinding:
      position: -1
  deas_dir:
    type: Directory
    inputBinding:
      position: 1
  genesets:
    type: File
    inputBinding:
      position: 2

