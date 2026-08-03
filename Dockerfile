# ECMExplainer — pinned runtime.
#
#   docker build -t ecmexplainer .
#   docker run --rm -p 3838:3838 ecmexplainer                  # local, no limits
#   docker run --rm -p 3838:3838 -e ECMX_MODE=web ecmexplainer # public instance
#
# rocker/r-ver pins both R and a dated CRAN snapshot, so the image rebuilds
# identically. Only five CRAN packages are needed — the deployed method requires
# no Bioconductor stack, because sign features need no batch correction.
FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev libpng-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
      shiny ggplot2 ggrepel commonmark \
    && rm -rf /tmp/downloaded_packages

WORKDIR /srv/ecmexplainer
COPY app.R METHODS.md ./
COPY R/    R/
COPY data/ data/

# Unprivileged: the app never writes to disk, so it needs no write access.
RUN useradd -r -u 10001 ecmx && chown -R ecmx:ecmx /srv/ecmexplainer
USER ecmx

EXPOSE 3838
ENV ECMX_MODE=web
HEALTHCHECK --interval=30s --timeout=8s --start-period=40s --retries=3 \
  CMD ["R", "-q", "-s", "-e", "q(status = if (inherits(try(readLines('http://127.0.0.1:3838', n = 1L, warn = FALSE), silent = TRUE), 'try-error')) 1L else 0L)"]

CMD ["R", "-q", "-e", "shiny::runApp('/srv/ecmexplainer', host='0.0.0.0', port=3838, launch.browser=FALSE)"]
