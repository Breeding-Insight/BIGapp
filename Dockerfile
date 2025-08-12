FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

# (Optional but recommended) freeze CRAN to "latest as of today" for reproducibility
# Change this date whenever you want to refresh to newer versions.
ENV CRAN_SNAPSHOT=2025-08-11

# System deps (compilers, headers, fonts; add wget for healthcheck)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran cmake git wget \
    libz-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libpng-dev libjpeg-dev libfreetype6-dev libfontconfig1-dev \
 && rm -rf /var/lib/apt/lists/*

# Make CRAN point to the snapshot while keeping Bioconductor repos
RUN R -q -e 'install.packages(c("BiocManager","remotes"), repos="https://cloud.r-project.org"); \
             local({ \
               r <- BiocManager::repositories(); \
               r["CRAN"] <- sprintf("https://packagemanager.posit.co/cran/%s", Sys.getenv("CRAN_SNAPSHOT","latest")); \
               options(repos = r); \
               cat("options(repos = ", deparse(r), ")\n", file="/usr/local/lib/R/etc/Rprofile.site", append=TRUE); \
             }); \
             BiocManager::install(version="3.20", ask=FALSE)'


# --- Bioconductor deps needed by BIGapp ---
RUN R -q -e 'BiocManager::install(c("Rsamtools","VariantAnnotation"), update=FALSE, ask=FALSE, Ncpus=2)'

# (Optional) speed up compiles for any source builds
ENV MAKEFLAGS="-j2"

# ----- Cache layer: install CRAN deps at pinned versions -----
WORKDIR /app
COPY DESCRIPTION /app/
# If you have NAMESPACE, keep this too:
# COPY NAMESPACE /app/

# ---- Resolve latest versions (at the snapshot) and install with install_version() ----
RUN Rscript - <<'RS'
options(Ncpus = 2)
pkgs <- c(
  "adegenet","curl","DT","dplyr","vcfR",
  "ggplot2","tidyr","shiny","config","bs4Dash",
  "golem","purrr","markdown","scales","plotly",
  "shinyWidgets","shinyjs","shinydisconnect","shinyalert",
  "stringr","updog","AGHmatrix","factoextra",
  "httr","future","shinycssloaders","RColorBrewer",
  "tibble","rrBLUP","MASS","Matrix","matrixcalc",
  "BIGr"
)

# Resolve exact versions from current repos (CRAN snapshot + Bioc)
ap <- available.packages(repos = getOption("repos"))
vers <- setNames(ap[pkgs, "Version"], pkgs)

# Print and save a copy-pasteable pin list
pins <- sprintf('%s@%s', names(vers), vers)
cat('Resolved pins:\n', paste(pins, collapse=',\n '), '\n')
dir.create('/opt/pins', recursive = TRUE, showWarnings = FALSE)
writeLines(pins, '/opt/pins/cran_pins.txt')

# Install each exact version (and skip vignettes/manual to speed up)
for (p in names(vers)) {
  remotes::install_version(
    p, version = vers[[p]], upgrade = "never",
    repos = getOption("repos"),
    INSTALL_opts = c("--no-build-vignettes","--no-manual")
  )
}
RS

# GitHub deps (pin to a commit if you can)
# ARG GITHUB_PAT
# ENV GITHUB_PAT=${GITHUB_PAT}
RUN Rscript -e 'remotes::install_github("jendelman/GWASpoly", upgrade="never", INSTALL_opts=c("--no-build-vignettes","--no-manual"))'


# ----- App source + install -----
COPY . /app
RUN R -q -e 'remotes::install_local(".", upgrade="never", INSTALL_opts=c("--no-build-vignettes","--no-manual"))'

# Security: run as non-root
RUN useradd -m appuser
USER appuser
WORKDIR /app

EXPOSE 80
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
  CMD wget -qO- http://localhost:${PORT:-80}/ || exit 1

# Read $PORT if provided; default to 80 (exec form recommended)
CMD ["R","-q","-e","options(shiny.port=as.integer(Sys.getenv('PORT','80')), shiny.host='0.0.0.0'); BIGapp::run_app()"]
