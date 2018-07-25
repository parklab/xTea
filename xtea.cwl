---
  hints: 
    - 
      dockerPull: "duplexa/xtea:v1"
      class: "DockerRequirement"
  arguments: []
  class: "CommandLineTool"
  inputs: 
    - 
      type: 
        - "File"
      id: "#bam"
      inputBinding: 
        position: 1
        prefix: "-b"
        separate: true
      secondaryFiles: 
        - "$(self.basename + '.bai')"
    - 
      type: 
        - "File"
      id: "#genome_tar"
      inputBinding: 
        position: 2
        prefix: "-l"
        separate: true
    - 
    - 
      type: 
        - "File"
      id: "#rep_libary_tar"
      inputBinding: 
        position: 3
        prefix: "-r"
        separate: true
    - 
      type: 
        - "int"
      id: "#nThreads"
      inputBinding: 
        position: 4
        prefix: "-n"
        separate: true
      default: 8
    - 
      type: 
        - "string"
      id: "#outdir"
      inputBinding: 
        position: 5
        prefix: "-p"
        separate: true
      default: "."
    - 
      type: 
        - "int"
      id: "#nclip"
      inputBinding: 
        position: 6
        prefix: "--nclip"
        separate: true
      default: 4
    - 
      type: 
        - "int"
      id: "#cr"
      inputBinding: 
        position: 7
        prefix: "--cr"
        separate: true
      default: 2
    - 
      type: 
        - "int"
      id: "#nd"
      inputBinding: 
        position: 8
        prefix: "--nd"
        separate: true
      default: 5
    - 
      type: 
        - "int"
      id: "#nfclip"
      inputBinding: 
        position: 9
        prefix: "--nfclip"
        separate: true
      default: 3
    - 
      type: 
        - "int"
      id: "#nfdisc"
      inputBinding: 
        position: 10
        prefix: "--nfdisc"
        separate: true
      default: 5
    - 
      type: 
        - "int"
      id: "#flklen"
      inputBinding: 
        position: 11
        prefix: "--flklen"
        separate: true
      default: 3000
    - 
      type: 
        - "int"
      id: "#f"
      inputBinding: 
        position: 12
        prefix: "--f"
        separate: true
      default: 19
  outputs: 
    - 
      type: 
        - "null"
        - "File"
      id: "#output"
      outputBinding: 
        glob: "candidate_disc_filtered_cns.txt"
  baseCommand: 
    - ["python", "gnrt_pipeline_cloud.py","-D", "-o run_jobs.sh", "-x /usr/local/bin/"]
  requirements: 
    - 
      class: "InlineJavascriptRequirement"
  cwlVersion: "draft-3"
