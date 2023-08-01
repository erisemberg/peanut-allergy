# A mutation in *Themis* contributes to peanut-induced oral anaphylaxis in CC027 mice 

This document describes how to reproduce heritability and QTL analyses of various peanut allergy-related phenotypes in a backcross between CC0027 (reacts orally to peanut) and C3H (does not react orally to peanut). 

Environment prep
----------------

This project uses a Docker container to produce an environment similar to that used in the original analysis (e.g. R v4.2.1 and R package versions available on August 1, 2023). In order to run this container you will need [Docker](https://docs.docker.com/get-docker/) installed. 

Build the docker container:

```
docker build . -t pnut 
```

Run the docker container, opening a terminal session within the container:

```
docker run -e PASSWORD=pw123 --rm -v $(pwd):/home/rstudio/work -p 8787:8787 -it gxd /bin/bash
```

Navigate to the working directory: 

```
cd home/rstudio/work 
```

Prepare R/qtl file 
-----------------------

```
Rscript rqtl-file-gen-pnutbc.R
```

Candidate gene analysis
-----------------------

Count total genes in each QTL region:

Count candidate genes in each QTL region:  

Summarize candidate gene data in table:  