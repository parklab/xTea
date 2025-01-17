---
  hints: 
    - 
      dockerPull: "warbler/xtea_denovo:v0.3.3"
      class: "DockerRequirement"
  arguments: []
  class: "CommandLineTool"
  inputs: 
    - 
      type: 
        - "string"
      id: "#id"
      inputBinding: 
        position: 1
        prefix: "-i"
        separate: true
    - 
      type: 
        - "File"
      id: "#bam"
      inputBinding: 
        position: 2
        prefix: "-b"
        separate: true
      secondaryFiles: 
        - $(self.basename.replace(/m$/,'i'))
    - 
      type: 
        - "File"
      id: "#pa_bam"
      inputBinding: 
        position: 3
        prefix: "--pa"
        separate: true
      secondaryFiles: 
        - $(self.basename.replace(/m$/,'i'))
    - 
      type: 
        - "File"
      id: "#ma_bam"
      inputBinding: 
        position: 4
        prefix: "--ma"
        separate: true
      secondaryFiles: 
        - $(self.basename.replace(/m$/,'i'))
    - 
      type: 
        - "File"
      id: "#genome_tar"
      inputBinding: 
        position: 5
        prefix: "-r"
        separate: true
    - 
      type: 
        - "File"
      id: "#rep_library_tar"
      inputBinding: 
        position: 6
        prefix: "-l"
        separate: true
    - 
      type: 
        - "int"
      id: "#nThreads"
      inputBinding: 
        position: 7
        prefix: "-n"
        separate: true
      default: 8
    - 
      type: 
        - "int"
      id: "#f"
      inputBinding: 
        position: 8
        prefix: "-f"
        separate: true
      default: 5907
    - 
      type: 
        - "int"
      id: "#y"
      inputBinding: 
        position: 9
        prefix: "-y"
        separate: true
      default: 7
  outputs: 
    - 
      type: 
        - "File"
      id: "#output"
      outputBinding: 
        glob: $(inputs.id)
  baseCommand: 
    - "python"
    - "/usr/local/bin/gnrt_pipeline_denovo_cloud_v38.py"
    - "-o"
    - "run_jobs.sh"
    - "--xtea"
    - "/usr/local/bin/"
    - "-p"
    - "."
  requirements: 
    - 
      class: "InlineJavascriptRequirement"
  cwlVersion: "v1.0"
