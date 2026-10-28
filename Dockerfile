# ECMExplainer — pinned runtime.
#
#   docker build -t ecmexplainer .
#   docker run --rm -p 3838:3838 ecmexplainer          # -> http://localhost:3838/ecmexplainer/
#
# Built on rocker/shiny, which bundles Shiny Server. Shiny Server serves the app
# under a path prefix (/ecmexplainer) with its asset URLs rewritten correctly;
# bare `runApp` emits absolute paths and breaks behind a sub-path proxy. The old
# /ecm/ path is kept as a static redirect so existing links still work.
#
# Shiny Server (open source) runs one R process for all visitors, so retraining
# runs in background R workers (future + promises) to keep the app responsive.
#
# rocker images point R at Posit Public Package Manager, so package installs
# pull prebuilt Ubuntu binaries instead of compiling from source.
FROM rocker/shiny:4.4.2

RUN install2.r --error --skipinstalled ggrepel commonmark future promises readxl \
    && rm -rf /tmp/downloaded_packages

# The default sample apps would otherwise be served alongside ours.
RUN rm -rf /srv/shiny-server/*

WORKDIR /srv/shiny-server/ecmexplainer
COPY app.R ABOUT.md METHODS.md ./
COPY R/            R/
COPY data/         data/
COPY supplementary/ supplementary/
COPY tests/        tests/
COPY www/          www/

# /ecm/ -> /ecmexplainer/ (a directory without app.R is served as static files).
RUN mkdir -p /srv/shiny-server/ecm && printf '%s\n' \
    '<!doctype html><meta charset="utf-8"><title>ECMExplainer</title>' \
    '<meta http-equiv="refresh" content="0; url=../ecmexplainer/">' \
    '<link rel="canonical" href="../ecmexplainer/">' \
    '<p>ECMExplainer has moved to <a href="../ecmexplainer/">/ecmexplainer/</a>.</p>' \
    > /srv/shiny-server/ecm/index.html

# Public instance limits. Override at run time for an unrestricted deployment.
ENV ECMX_MODE=web \
    ECMX_MAX_UPLOAD_MB=50 \
    ECMX_MAX_SAMPLES=500 \
    ECMX_MAX_TRAIN=200 \
    ECMX_MAX_PER_PHENO=60 \
    ECMX_WORKERS=2 \
    ECMX_CORES=2

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && chown -R shiny:shiny /srv/shiny-server
USER shiny

EXPOSE 3838
HEALTHCHECK --interval=30s --timeout=8s --start-period=60s --retries=3 \
  CMD ["R", "-q", "-s", "-e", "q(status = if (inherits(try(readLines('http://127.0.0.1:3838/ecmexplainer/', n = 1L, warn = FALSE), silent = TRUE), 'try-error')) 1L else 0L)"]

CMD ["/usr/local/bin/docker-entrypoint.sh"]
