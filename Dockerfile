# syntax=docker/dockerfile:1.7
# Buildx/Actions will pass BASE_IMAGE as a manifest tag that covers both arches
ARG BASE_IMAGE=docker.io/breedinginsight/bigapp-deps:latest
FROM ${BASE_IMAGE}

SHELL ["/bin/bash","-eo","pipefail","-c"]
ENV DEBIAN_FRONTEND=noninteractive TZ=UTC
ENV MAKEFLAGS="-j2" R_PKG_INSTALL_ARGS="--no-build-vignettes --no-manual"

# App install (only your code changes should rebuild this layer)
WORKDIR /app
COPY DESCRIPTION /app/
# COPY NAMESPACE /app/  # if present, include for better cache hits
COPY . /app
RUN R -q -e "remotes::install_local('.', upgrade='never', dependencies=TRUE, \
                                    INSTALL_opts=c('--no-build-vignettes','--no-manual')); \
             library(BIGapp)" \
    || (echo 'ERROR: BIGapp or one of its dependencies failed to install' >&2; exit 1)

# Runtime
RUN useradd -m appuser
USER appuser
WORKDIR /app
EXPOSE 80
HEALTHCHECK --interval=30s --timeout=5s --retries=5 \
  CMD wget -qO- http://localhost:${PORT:-80}/ || exit 1
CMD ["R","-q","-e","options(shiny.port=as.integer(Sys.getenv('PORT','80')), shiny.host='0.0.0.0'); BIGapp::run_app()"]
