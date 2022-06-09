#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.gencode)
        entryname: gencode.annotation.gff3

hints:
  - class: DockerRequirement
    dockerPull: warbler/xtea_germline:v0.2.1

baseCommand:
  - "python"
  - "/usr/local/bin/gnrt_pipeline_germline_cloud_v38.py"
  - "-D"
  - "-o"
  - "run_jobs.sh"
  - "--xtea"
  - "/usr/local/bin/"
  - "-p"
  - "."

inputs:
  - id: sample
    type: string
    inputBinding:
      prefix: -i
    doc: sample name

  - id: input_bam
    type: File
    inputBinding:
      prefix: -b
    secondaryFiles:
      - .bai

  - id: genome_tar
    type: File
    inputBinding:
      prefix: -r
    doc: genome library |
         gzip archive

  - id: repeats_tar
    type: File
    inputBinding:
      prefix: -l
    doc: repeats library |
         gzip archive

  - id: gencode
    type: File
    doc: gencode.annotation.gff3

  - id: nthreads
    type: int
    default: 16
    inputBinding:
      prefix: -n

  - id: f
    type: int
    default: 5907
    inputBinding:
      prefix: -f

  - id: y
    type: int
    default: 7
    inputBinding:
      prefix: -y
    doc: analysis type |
         1 for LINE1, 2 for ALU, 4 for SVA

outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.sample + "_*" + ".tar.gz")

doc: |
  run xTea
