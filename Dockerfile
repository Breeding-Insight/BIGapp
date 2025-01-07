FROM rocker/r-ver:4.4.2
RUN apt-get update && apt-get install -y  cmake libz-dev libcurl4-openssl-dev libssl-dev
RUN R -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("adegenet",upgrade="never", version = "2.1.10")'
RUN Rscript -e 'remotes::install_version("curl",upgrade="never", version = "6.0.1")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.33")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.1.4")'
RUN Rscript -e 'remotes::install_version("vcfR",upgrade="never", version = "1.15.0")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.5.1")'
RUN Rscript -e 'remotes::install_version("tidyr",upgrade="never", version = "1.3.1")'
RUN Rscript -e 'remotes::install_version("curl",upgrade="never", version = "6.0.1")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.9.1")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.2")'
RUN Rscript -e 'remotes::install_version("bs4Dash",upgrade="never", version = "2.3.4")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.5.1")'
RUN Rscript -e 'remotes::install_version("purrr",upgrade="never", version = "1.0.2")'
RUN Rscript -e 'remotes::install_version("markdown",upgrade="never", version = "1.13")'
RUN Rscript -e 'remotes::install_version("scales",upgrade="never", version = "1.3.0")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.4")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.8.7")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.1.0")'
RUN Rscript -e 'remotes::install_version("shinydisconnect",upgrade="never", version = "0.1.1")'
RUN Rscript -e 'remotes::install_version("shinyalert",upgrade="never", version = "3.1.0")'
RUN Rscript -e 'remotes::install_version("stringr",upgrade="never", version = "1.5.1")'
RUN Rscript -e 'remotes::install_version("updog",upgrade="never", version = "2.1.5")'
RUN Rscript -e 'remotes::install_version("AGHmatrix",upgrade="never", version = "2.1.4")'
RUN Rscript -e 'remotes::install_version("factoextra",upgrade="never", version = "1.0.7")'
RUN Rscript -e 'remotes::install_version("httr",upgrade="never", version = "1.4.7")'
RUN Rscript -e 'remotes::install_version("future",upgrade="never", version = "1.34.0")'
RUN Rscript -e 'remotes::install_version("shinycssloaders",upgrade="never", version = "1.1.0")'
RUN Rscript -e 'remotes::install_version("RColorBrewer",upgrade="never", version = "1.1.3")'
RUN Rscript -e 'remotes::install_version("tibble",upgrade="never", version = "3.2.1")'
RUN Rscript -e 'remotes::install_version("rrBLUP",upgrade="never", version = "4.6.3")'
RUN Rscript -e 'remotes::install_version("MASS",upgrade="never", version = "7.3.60.2")'
RUN Rscript -e 'remotes::install_version("Matrix",upgrade="never", version = "1.7.0")'
RUN Rscript -e 'remotes::install_version("matrixcalc",upgrade="never", version = "1.0.6")'
RUN Rscript -e 'remotes::install_github("Breeding-Insight/BIGr",upgrade="never")'
RUN Rscript -e 'remotes::install_github("jendelman/GWASpoly",upgrade="never")'

RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /build_zone
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');BIGapp::run_app()"

