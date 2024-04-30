FROM rocker/shiny:latest

# Install build dependencies, liblzma-dev and libbz2-dev for Rhtslib and libglpk40 for adegenet.
RUN sudo apt update
RUN sudo apt -y install build-essential liblzma-dev libbz2-dev libglpk40

# Use all package repositories. Install app dependencies.
RUN R -e 'setRepositories(ind=1:7); install.packages(c("updog", "ggplot2", "VariantAnnotation", "SNPRelate", "adegenet", "future", "scales", "AGHmatrix", "stats", "factoextra", "readxl", "ggrepel", "dplyr", "shiny", "shinydashboard","randomcoloR","plotly", "DT","RColorBrewer", "dichromat", "bs4Dash", "shinyWidgets","data.table", "matrixcalc","Matrix", "shinyalert","rrBLUP", "tidyverse"))'

# "GWASpoly" is not in any of the repositories, install from GitHub (https://github.com/jendelman/GWASpoly/).
RUN R -e 'install.packages("devtools"); devtools::install_github("jendelman/GWASpoly", build_vignettes=FALSE)'

# Clean up the shiny server example apps.
RUN sudo rm -rf /srv/shiny-server/

# Copy app source code.
COPY  ./BIG_app/ /srv/shiny-server/app/

# Copy shiny server config.
COPY ./shiny-server.conf /etc/shiny-server/shiny-server.conf

CMD ["/usr/bin/shiny-server"]