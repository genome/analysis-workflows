[![Build Status](https://travis-ci.org/genome/analysis-workflows.svg?branch=master)](https://travis-ci.org/genome/analysis-workflows)

# analysis-workflows

## Overview
The [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) and [contributing](https://github.com/genome/analysis-workflows#contributions) staff, faculty, labs and departments of [Washington University School of Medicine](https://medicine.wustl.edu/) (WUSM) share [Common Workflow Language](https://www.commonwl.org/) (CWL) workflow defintions focused on reusable, reproducible analysis pipelines for genomics data.  

## Structure

The main structure of this repo is described in the following table:

| Path | Description |
| --- | --- |
| definitions | parent directory containing all CWL tool and workflow definitions |
| definitions/pipelines | all workflows which rely on subworkflows and tools to produce final outputs |
| definitions/subworkflows | workflows that combine multiple tools to produce intermediate (used as inputs to other subworkflows) pipeline outputs |
| definitions/tools | CWL that wrap command line interfaces or scripts connecting multiple tools |
| definitions/types | custom CWL data types for inputs to tools and workflows |
| example_data | example input data, input YAML files, and expected output files for testing |

## Documentation

All documentation of CWL pipelines, subworkflows, and tools as well as additional information regarding test data, continous integration, and configuration can be found on the GitHub wiki:
https://github.com/genome/analysis-workflows/wiki

## Images

All MGI supported Docker images used in the tool workflow definitions are available on [mgibio DockerHub](https://hub.docker.com/u/mgibio/): 

Many tools rely on third-party Docker images publicly available from sources such as [Docker Hub](https://hub.docker.com) and [BioContainers](https://biocontainers.pro).

## Installation

No software installation is necessary to use these CWL definitions. 

Each CWL file is validated using [cwltool](https://github.com/common-workflow-language/cwltool). Additionl workflow definition testing is performed with [Cromwell](https://github.com/broadinstitute/cromwell). However, currently there are no automated workflow tests using Cromwell.

### Workflow Execution Service

These workflow definitions are built for interoperability with any [Workflow Execution Service](https://github.com/ga4gh/workflow-execution-service-schemas) (WES) schema compatible implementation that supports CWL.

## Contributions

A big thanks to all of the developers, bioinformaticians and scientists who built this resource. For a complete list of software contributions, i.e. commits, to this repository, please see the [GitHub Contributors](https://github.com/genome/analysis-workflows/graphs/contributors) section.

### Collaborators

The following WUSM collaborators have provided significant contributions in terms of workflow design, scientific direction, and validation of [analysis-workflows](https://github.com/genome/analysis-workflows) output.

* [Chris Miller, PhD](https://www.genome.wustl.edu/people/chris-miller-phd/)
* [David Spencer, MD, PhD](https://www.genome.wustl.edu/people/david-spencer/)
* [Malachi Griffith, PhD](https://www.genome.wustl.edu/people/malachi-griffith/)
* [Obi L. Griffith, PhD](https://www.genome.wustl.edu/people/obi-griffith/)

### Departments, Institutes, and Labs
* [Department of Genetics](http://genetics.wustl.edu/)
* [Department of Medicine, Divisions of Hematology and Oncology](https://oncology.wustl.edu/)
* [Griffith Lab](https://www.genome.wustl.edu/research/labs/griffith-lab/)
* [Institue of Clinical and Translational Sciences](https://icts.wustl.edu/)
* [McDonnell Genome Institute](https://www.genome.wustl.edu/)


[![DOI](https://zenodo.org/badge/64162512.svg)](https://zenodo.org/badge/latestdoi/64162512)

