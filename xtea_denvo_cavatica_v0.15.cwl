{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "baseCommand": [
        "python",
        "/usr/local/bin/gnrt_pipeline_denovo_cloud_v38.py",
        "-o",
        "run_jobs.sh",
        "--xtea",
        "/usr/local/bin/",
        "-p",
        "."
    ],
    "inputs": [
        {
            "id": "#id",
            "type": "string",
            "inputBinding": {
                "prefix": "-i",
                "shellQuote": false,
                "position": 1
            }
        },
        {
            "id": "#bam",
            "type": "File",
            "inputBinding": {
                "prefix": "-b",
                "shellQuote": false,
                "position": 2
            },
            "secondaryFiles": [
                {
                    "pattern": "$(self.basename.replace(/m$/,'i'))",
                    "required": true
                }
            ]
        },
        {
            "id": "#pa_bam",
            "type": "File",
            "inputBinding": {
                "prefix": "--pa",
                "shellQuote": false,
                "position": 3
            },
            "secondaryFiles": [
                {
                    "pattern": "$(self.basename.replace(/m$/,'i'))",
                    "required": true
                }
            ]
        },
        {
            "id": "#ma_bam",
            "type": "File",
            "inputBinding": {
                "prefix": "--ma",
                "shellQuote": false,
                "position": 4
            },
            "secondaryFiles": [
                {
                    "pattern": "$(self.basename.replace(/m$/,'i'))",
                    "required": true
                }
            ]
        },
        {
            "id": "#genome_tar",
            "type": "File",
            "inputBinding": {
                "prefix": "-r",
                "shellQuote": false,
                "position": 5
            }
        },
        {
            "id": "#rep_library_tar",
            "type": "File",
            "inputBinding": {
                "prefix": "-l",
                "shellQuote": false,
                "position": 6
            }
        },
        {
            "id": "#nThreads",
            "type": "int",
            "inputBinding": {
                "prefix": "-n",
                "shellQuote": false,
                "position": 7
            },
            "default": 8
        },
        {
            "id": "#f",
            "type": "int",
            "inputBinding": {
                "prefix": "-f",
                "shellQuote": false,
                "position": 8
            },
            "default": 5907
        },
        {
            "id": "#y",
            "type": "int",
            "inputBinding": {
                "prefix": "-y",
                "shellQuote": false,
                "position": 9
            },
            "default": 7
        },
        {
            "id": "#decompress",
            "type": "string",
            "inputBinding": {
                "prefix": "-D",
                "shellQuote": false,
                "position": 10
            },
            "default": " "
        }
    ],
    "outputs": [
        {
            "id": "#output",
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.id)"
            }
        }
    ],
    "label": "xtea-denovo-cwl",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": 24,
            "coresMin": 0
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "hints": [
        {
            "class": "DockerRequirement",
            "dockerPull": "warbler/xtea_denovo:v0.4.3"
        }
    ],
    "sbg:projectName": "KidsFirst_denovo_MEI",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1637161775,
            "sbg:revisionNotes": "Copy of simon_chu/test-workflow/xtea-denovo-cwl/13"
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1637164414,
            "sbg:revisionNotes": "v0.1.9"
        },
        {
            "sbg:revision": 2,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646540374,
            "sbg:revisionNotes": "v0.3.1"
        },
        {
            "sbg:revision": 3,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646572005,
            "sbg:revisionNotes": "revise_output"
        },
        {
            "sbg:revision": 4,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646572434,
            "sbg:revisionNotes": "test4"
        },
        {
            "sbg:revision": 5,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646572571,
            "sbg:revisionNotes": "test5"
        },
        {
            "sbg:revision": 6,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646596924,
            "sbg:revisionNotes": "test7"
        },
        {
            "sbg:revision": 7,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646752959,
            "sbg:revisionNotes": "update purity"
        },
        {
            "sbg:revision": 8,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646774757,
            "sbg:revisionNotes": "with_purity"
        },
        {
            "sbg:revision": 9,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646793821,
            "sbg:revisionNotes": "debug version"
        },
        {
            "sbg:revision": 10,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646801527,
            "sbg:revisionNotes": "beta"
        },
        {
            "sbg:revision": 11,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646856039,
            "sbg:revisionNotes": "test_output"
        },
        {
            "sbg:revision": 12,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1646863839,
            "sbg:revisionNotes": "latest"
        },
        {
            "sbg:revision": 13,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1650222016,
            "sbg:revisionNotes": "pseudogene"
        },
        {
            "sbg:revision": 14,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1650229818,
            "sbg:revisionNotes": "pseudogene2"
        },
        {
            "sbg:revision": 15,
            "sbg:modifiedBy": "simon_chu",
            "sbg:modifiedOn": 1650230964,
            "sbg:revisionNotes": "Pseudogene3"
        }
    ],
    "sbg:image_url": null,
    "sbg:appVersion": [
        "v1.2"
    ],
    "id": "https://cavatica-api.sbgenomics.com/v2/apps/simon_chu/kidsfirst-denovo-mei/xtea-denovo-cwl/15/raw/",
    "sbg:id": "simon_chu/kidsfirst-denovo-mei/xtea-denovo-cwl/15",
    "sbg:revision": 15,
    "sbg:revisionNotes": "Pseudogene3",
    "sbg:modifiedOn": 1650230964,
    "sbg:modifiedBy": "simon_chu",
    "sbg:createdOn": 1637161775,
    "sbg:createdBy": "simon_chu",
    "sbg:project": "simon_chu/kidsfirst-denovo-mei",
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [],
    "sbg:contributors": [
        "simon_chu"
    ],
    "sbg:latestRevision": 15,
    "sbg:publisher": "sbg",
    "sbg:content_hash": "a3ac3d347c6d9d417f9fc243da7877599c96058b60f22d208d1f6574a9e64694f",
    "sbg:workflowLanguage": "CWL"
}