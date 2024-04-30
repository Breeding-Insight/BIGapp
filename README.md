# Genomics_Shiny_App

Creating a shiny app that includes the mature genomics/bioinformatics pipelines

To use:
1. install.packages("shiny")
2. library("shiny")
3. setwd("/Genomics_Shiny_App")
4. runApp("BIG_app")
5. View shiny app in browser


# Docker
## Build docker image
> This may take 20 minutes or so to build the first time, as many dependencies need to be downloaded and built.

> Do not neglect the trailing `.`, this tells docker to look for a Dockerfile in the current working directory.

On x86 (e.g. intel) use:
`docker build -t shinyimage .`

On ARM (e.g. apple silicon macs), use:
`DOCKER_DEFAULT_PLATFORM=linux/amd64 docker build -t shinyimage .`


## Run docker image
`docker run -p 3838:3838 shinyimage`

## Access in your browser
http://localhost:3838/