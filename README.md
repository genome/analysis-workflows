# cancer-genomics-workflow
CWL workflow defintions for executing MGI Cancer Genomics analysis.

## Docker Images

At this time, the toil_compatibility branch relies on large, aka. fat, Docker images. Below is a mapping of workflow to Docker image name, all available on [DockerHub](https://hub.docker.com/u/mgibio/)

| Workflow | Docker Image | Main CWL |
| --- | --- | --- |
| unaligned_bam_to_bqsr | mgibio/dna-alignment | unaligned_bam_to_bgsr/workflow.cwl |
| detect_variants | mgibio/cle | detect_variants/detect_variants.cwl |
| rnaseq | mgibio/rnaseq | rnaseq/workflow.cwl |
| tumor_only_exome | mgibio/cle | exome_workflow.cwl |

If a tool or workflow is a subworkflow of a larger workflow, then the same Docker image is required.

