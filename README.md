# analysis-workflows

This repo contains CWL workflow defintions for executing MGI analysis pipelines. The structure of the repo is as follows:
| Path | Description |
| --- | --- |
| definitions | parent directory containing all CWL tool and workflow definitions |
| definitions/pipelines | all workflows which rely on subworkflows and tools to produce final outputs |
| definitions/subworkflows | workflows that combine multiple tools to produce intermediate (used as inputs to other subworkflows) pipeline outputs |
| definitions/tools | CWL that wrap command line interfaces or scripts connecting multiple tools |
| example_data | Example input data, inputs YAML files, and expected output files |

## Docker Images

The master branch specifies `DockerRequirement` for each tool in the workflow. All available on [DockerHub](https://hub.docker.com/u/mgibio/):

| Workflow | Docker Image | Main CWL |
| --- | --- | --- |
| unaligned_bam_to_bqsr | mgibio/dna-alignment | unaligned_bam_to_bqsr/workflow.cwl |
| detect_variants | mgibio/cle | detect_variants/detect_variants.cwl |
| rnaseq | mgibio/rnaseq | rnaseq/workflow.cwl |
| tumor_only_exome | mgibio/cle | exome_workflow.cwl |
| umi_alignment | mgibio/dna-alignment | umi_alignment/duplex_workflow.cwl |

If a tool or workflow is a subworkflow of a larger workflow, then the same Docker image is required.



[![DOI](https://zenodo.org/badge/64162512.svg)](https://zenodo.org/badge/latestdoi/64162512)

