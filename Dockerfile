# syntax=docker/dockerfile:1.7
FROM rocker/r2u:noble

SHELL ["/bin/bash","-eo","pipefail","-c"]
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# keep compiles tame on GH runners / qemu
ENV MAKEFLAGS="-j2"
ENV R_PKG_INSTALL_ARGS="--no-build-vignettes --no-manual"

# --- System deps (no Bioc via apt) ---
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran cmake git wget pkg-config \
      ccache \
      libhts-dev \
      libcurl4-gnutls-dev libssl-dev libxml2-dev \
      libpng-dev libjpeg-dev libfreetype6-dev libfontconfig1-dev \
      libbz2-dev liblzma-dev zlib1g-dev \
      r-cran-remotes r-cran-pak r-cran-biocmanager \
    && rm -rf /var/lib/apt/lists/*

# --- ccache for repeated builds ---
ENV CCACHE_DIR=/root/.cache/ccache
RUN --mount=type=cache,target=/root/.cache/ccache \
    ccache -M 2G && \
    { echo 'CC=ccache gcc'; echo 'CXX=ccache g++'; echo 'FC=ccache gfortran'; } >> /usr/lib/R/etc/Makevars.site

# --- CRAN via r2u binaries first (fast) ---
ENV CRAN_PKGS="adegenet curl DT dplyr vcfR \
ggplot2 tidyr shiny config bs4Dash \
golem purrr markdown scales plotly \
shinyWidgets shinyjs shinydisconnect shinyalert \
stringr updog AGHmatrix factoextra \
httr future shinycssloaders RColorBrewer \
tibble rrBLUP MASS Matrix matrixcalc BIGr"

RUN install2.r --skipinstalled --ncpus 1 $CRAN_PKGS || true

# --- Bioconductor from source (sequential; correct for R 4.5 -> Bioc 3.21) ---
RUN --mount=type=cache,target=/root/.cache/R/src Rscript - <<'RS'
options(Ncpus = 2)
repos <- BiocManager::repositories(version = "3.21")
options(repos = repos)
bioc_pkgs <- c("Rsamtools","Biostrings","pwalign","VariantAnnotation")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing Bioconductor package from source: ", p)
    install.packages(p,
      dependencies = TRUE,
      type = "source",
      INSTALL_opts = c("--no-build-vignettes","--no-manual")
    )
  }
}
RS

# --- Fallback for any missing CRAN pkgs (sequential, no pak helper) ---
RUN --mount=type=cache,target=/root/.cache/R/src Rscript - <<'RS'
pkgs <- strsplit(Sys.getenv("CRAN_PKGS"), " +")[[1]]
miss <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss)) {
  message("Installing remaining CRAN from source: ", paste(miss, collapse=", "))
  for (p in miss) {
    install.packages(p,
      repos = "https://cloud.r-project.org",
      type  = "source",
      INSTALL_opts = c("--no-build-vignettes","--no-manual")
    )
  }
}
miss2 <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss2)) {
  message("Missing after fallback: ", paste(miss2, collapse=", "))
  quit(status = 1)
}
RS

# --- GitHub dep (avoid pak helper under QEMU) ---
# ARG GITHUB_PAT
# ENV GITHUB_PAT=${GITHUB_PAT}
RUN R -q -e "remotes::install_github('jendelman/GWASpoly', upgrade='never', dependencies=TRUE, \
                                     INSTALL_opts=c('--no-build-vignettes','--no-manual'))"

# --- App install (cache-friendly) ---
WORKDIR /app
COPY DESCRIPTION /app/
# COPY NAMESPACE /app/   # uncomment if present

RUN R -q -e "pkgs <- c('adegenet','curl','DT','dplyr','vcfR','ggplot2','tidyr','shiny','config','bs4Dash', \
                       'golem','purrr','markdown','scales','plotly','shinyWidgets','shinyjs','shinydisconnect','shinyalert', \
                       'stringr','updog','AGHmatrix','factoextra','httr','future','shinycssloaders','RColorBrewer', \
                       'tibble','rrBLUP','MASS','Matrix','matrixcalc','BIGr', \
                       'Rsamtools','VariantAnnotation','Biostrings','pwalign'); \
                 ip <- installed.packages(); dir.create('/opt/pins', recursive=TRUE, showWarnings=FALSE); \
                 write.table(data.frame(Package=pkgs, Version=ip[pkgs,'Version']), \
                             file='/opt/pins/installed_versions.tsv', sep='\t', row.names=FALSE, quote=FALSE)"

COPY . /app
RUN --mount=type=cache,target=/root/.cache/R/src \
    R -q -e "remotes::install_local('.', upgrade='never', dependencies=TRUE, \
                                    INSTALL_opts=c('--no-build-vignettes','--no-manual'))"

# --- Runtime ---
RUN useradd -m appuser
USER appuser
WORKDIR /app

EXPOSE 80
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
  CMD wget -qO- http://localhost:${PORT:-80}/ || exit 1

CMD ["R","-q","-e","options(shiny.port=as.integer(Sys.getenv('PORT','80')), shiny.host='0.0.0.0'); BIGapp::run_app()"]
