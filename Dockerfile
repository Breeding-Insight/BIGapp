# syntax=docker/dockerfile:1.7
FROM rocker/r2u:noble

SHELL ["/bin/bash","-eo","pipefail","-c"]
ENV DEBIAN_FRONTEND=noninteractive TZ=UTC
# Keep compiles modest to avoid OOM on small runners
ENV MAKEFLAGS="-j2"
# Skip vignettes/manuals for speed
ENV R_PKG_INSTALL_ARGS="--no-build-vignettes --no-manual"

# ---------- System deps ----------
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran cmake git wget \
      ccache \
      libhts-dev \
      libcurl4-gnutls-dev libssl-dev libxml2-dev \
      libpng-dev libjpeg-dev libfreetype6-dev libfontconfig1-dev \
      libbz2-dev liblzma-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------- ccache ----------
ENV CCACHE_DIR=/root/.cache/ccache
RUN --mount=type=cache,target=/root/.cache/ccache \
    ccache -M 2G && \
    { echo 'CC=ccache gcc'; echo 'CXX=ccache g++'; echo 'FC=ccache gfortran'; } >> /usr/lib/R/etc/Makevars.site

# ---------- pak + BiocManager (prefer apt, fall back to CRAN) ----------
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends r-cran-pak r-cran-biocmanager || true
RUN R -q -e 'if (!"pak" %in% rownames(installed.packages())) install.packages("pak", repos="https://cloud.r-project.org")'
RUN R -q -e 'if (!"BiocManager" %in% rownames(installed.packages())) install.packages("BiocManager", repos="https://cloud.r-project.org")'

# ---------- CRAN via binaries first ----------
ENV CRAN_PKGS="adegenet curl DT dplyr vcfR \
ggplot2 tidyr shiny config bs4Dash \
golem purrr markdown scales plotly \
shinyWidgets shinyjs shinydisconnect shinyalert \
stringr updog AGHmatrix factoextra \
httr future shinycssloaders RColorBrewer \
tibble rrBLUP MASS Matrix matrixcalc BIGr"

# Try r2u binaries (no --error, low parallelism)
RUN install2.r --skipinstalled --ncpus 1 $CRAN_PKGS || true

# ---------- Fallback + final check in one R call (no quoting issues) ----------
RUN --mount=type=cache,target=/root/.cache/R/pak Rscript - <<'RS'
pkgs <- strsplit(Sys.getenv("CRAN_PKGS"), " +")[[1]]
installed <- rownames(installed.packages())
miss <- setdiff(pkgs, installed)

if (length(miss)) {
  message("Installing from source with pak, sequentially: ", paste(miss, collapse=", "))
  for (p in miss) {
    pak::pkg_install(p, ask = FALSE)
  }
}
installed2 <- rownames(installed.packages())
miss2 <- setdiff(pkgs, installed2)
if (length(miss2)) {
  message("Missing after fallback: ", paste(miss2, collapse=", "))
  quit(status = 1)
}
RS

# ---------- Bioconductor (compiled; libhts-dev helps) ----------
RUN --mount=type=cache,target=/root/.cache/R/pak \
    R -q -e "pak::repo_add(bioc = BiocManager::repositories()[['BioCsoft']]); \
             pak::pkg_install(c('bioc::Rsamtools','bioc::VariantAnnotation'), ask=FALSE)"

# ---------- GitHub dep (pin to a commit if you can) ----------
# ARG GITHUB_PAT
# ENV GITHUB_PAT=${GITHUB_PAT}
RUN --mount=type=cache,target=/root/.cache/R/pak \
    R -q -e "pak::pkg_install('github::jendelman/GWASpoly', ask=FALSE)"

# ---------- App (cache-friendly) ----------
WORKDIR /app
COPY DESCRIPTION /app/
# COPY NAMESPACE /app/   # uncomment if present for better cache hits

# Snapshot versions (optional but handy)
RUN R -q -e "pkgs <- c('adegenet','curl','DT','dplyr','vcfR','ggplot2','tidyr','shiny','config','bs4Dash', \
                       'golem','purrr','markdown','scales','plotly','shinyWidgets','shinyjs','shinydisconnect','shinyalert', \
                       'stringr','updog','AGHmatrix','factoextra','httr','future','shinycssloaders','RColorBrewer', \
                       'tibble','rrBLUP','MASS','Matrix','matrixcalc','BIGr', \
                       'Rsamtools','VariantAnnotation'); \
                 ip <- installed.packages(); dir.create('/opt/pins', recursive=TRUE, showWarnings=FALSE); \
                 write.table(data.frame(Package=pkgs, Version=ip[pkgs,'Version']), \
                             file='/opt/pins/installed_versions.tsv', sep='\t', row.names=FALSE, quote=FALSE)"

COPY . /app
RUN --mount=type=cache,target=/root/.cache/R/pak \
    R -q -e "pak::pkg_install('local::.', ask=FALSE)"

# ---------- Runtime ----------
RUN useradd -m appuser
USER appuser
WORKDIR /app

EXPOSE 80
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
  CMD wget -qO- http://localhost:${PORT:-80}/ || exit 1

CMD ["R","-q","-e","options(shiny.port=as.integer(Sys.getenv('PORT','80')), shiny.host='0.0.0.0'); BIGapp::run_app()"]
