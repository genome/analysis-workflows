# analysis-workflows

This repo contains CWL workflow defintions for executing MGI analysis pipelines. 

The structure of the repo is as follows:

| Path | Description |
| --- | --- |
| definitions | parent directory containing all CWL tool and workflow definitions |
| definitions/pipelines | all workflows which rely on subworkflows and tools to produce final outputs |
| definitions/subworkflows | workflows that combine multiple tools to produce intermediate (used as inputs to other subworkflows) pipeline outputs |
| definitions/tools | CWL that wrap command line interfaces or scripts connecting multiple tools |
| definitions/types | custom CWL data types for inputs to tools and workflows |
| example_data | example input data, input YAML files, and expected output files for testing |

## CWL Documentation

For all documentation of pipelines, subworkflows and tools as well as additional information regarding test data, continous integration and configuration, please see the GitHub wiki:
https://github.com/genome/analysis-workflows/wiki

## Docker Images

All MGI supported Docker images used in the tool workflow definitions are available on [mgibio DockerHub](https://hub.docker.com/u/mgibio/): 

Many tools rely on third-party Docker images publicly available from sources such as [Docker Hub](https://hub.docker.com) and [BioContainers](https://biocontainers.pro).


[![DOI](https://zenodo.org/badge/64162512.svg)](https://zenodo.org/badge/latestdoi/64162512)

