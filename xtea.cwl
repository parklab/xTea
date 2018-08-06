---
  hints: 
    - 
      dockerPull: "duplexa/xtea:v2"
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
        - $(self.basename.replace(/m$/,'i'))
    - 
      type: 
        - "File"
      id: "#genome_tar"
      inputBinding: 
        position: 2
        prefix: "-l"
        separate: true
    - 
      type: 
        - "File"
      id: "#rep_library_tar"
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
        - "int"
      id: "#nclip"
      inputBinding: 
        position: 5
        prefix: "--nclip"
        separate: true
      default: 4
    - 
      type: 
        - "int"
      id: "#cr"
      inputBinding: 
        position: 6
        prefix: "--cr"
        separate: true
      default: 2
    - 
      type: 
        - "int"
      id: "#nd"
      inputBinding: 
        position: 7
        prefix: "--nd"
        separate: true
      default: 5
    - 
      type: 
        - "int"
      id: "#nfclip"
      inputBinding: 
        position: 8
        prefix: "--nfclip"
        separate: true
      default: 3
    - 
      type: 
        - "int"
      id: "#nfdisc"
      inputBinding: 
        position: 9
        prefix: "--nfdisc"
        separate: true
      default: 5
    - 
      type: 
        - "int"
      id: "#flklen"
      inputBinding: 
        position: 10
        prefix: "--flklen"
        separate: true
      default: 3000
    - 
      type: 
        - "int"
      id: "#f"
      inputBinding: 
        position: 11
        prefix: "-f"
        separate: true
      default: 19
  outputs: 
    - 
      type: 
        - "File"
      id: "#output"
      outputBinding: 
        glob: "results/results.tar.gz"
  baseCommand: 
    - "python"
    - "/usr/local/bin/gnrt_pipeline_cloud.py"
    - "-D"
    - "-o"
    - "run_jobs.sh"
    - "-x"
    - "/usr/local/bin/"
    - "-p"
    - "."
  requirements: 
    - 
      class: "InlineJavascriptRequirement"
  cwlVersion: "v1.0"
