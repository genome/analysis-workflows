# analysis-workflows

CWL workflow defintions for executing MGI analysis pipelines.

## Docker Images

The master branch specifies `DockerRequirement` for each tool in the workflow.

By contrast, the toil_compatibility branch relies on large, a/k/a fat, Docker images. Below is a mapping of workflow to Docker image name, all available on [DockerHub](https://hub.docker.com/u/mgibio/):

| Workflow | Docker Image | Main CWL |
| --- | --- | --- |
| unaligned_bam_to_bqsr | mgibio/dna-alignment | unaligned_bam_to_bqsr/workflow.cwl |
| detect_variants | mgibio/cle | detect_variants/detect_variants.cwl |
| rnaseq | mgibio/rnaseq | rnaseq/workflow.cwl |
| tumor_only_exome | mgibio/cle | exome_workflow.cwl |
| umi_alignment | mgibio/dna-alignment | umi_alignment/duplex_workflow.cwl |

If a tool or workflow is a subworkflow of a larger workflow, then the same Docker image is required.



[![DOI](https://zenodo.org/badge/64162512.svg)](https://zenodo.org/badge/latestdoi/64162512)

