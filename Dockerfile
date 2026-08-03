# ECMExplainer — pinned runtime.
#
#   docker build -t ecmexplainer .
#   docker run --rm -p 3838:3838 ecmexplainer          # -> http://localhost:3838/ecm/
#
# Built on rocker/shiny, which bundles Shiny Server. Two reasons:
#   1. Shiny Server serves the app under a path prefix (/ecm) with its asset URLs
#      rewritten correctly; bare `runApp` emits absolute /shared/... paths and
#      breaks behind a sub-path proxy.
#   2. It spawns one R process per session, so a slow retrain by one visitor
#      does not block everyone else.
#
# rocker images point R at Posit Public Package Manager, so package installs
# pull prebuilt Ubuntu binaries instead of compiling from source.
FROM rocker/shiny:4.4.2

RUN install2.r --error --skipinstalled ggrepel commonmark \
    && rm -rf /tmp/downloaded_packages

# The default sample apps would otherwise be served alongside ours.
RUN rm -rf /srv/shiny-server/*

WORKDIR /srv/shiny-server/ecm
COPY app.R METHODS.md ./
COPY R/    R/
COPY data/ data/

# Public instance limits. Override at run time for an unrestricted deployment.
ENV ECMX_MODE=web \
    ECMX_MAX_UPLOAD_MB=20 \
    ECMX_MAX_SAMPLES=60 \
    ECMX_MAX_TRAIN=60

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && chown -R shiny:shiny /srv/shiny-server
USER shiny

EXPOSE 3838
HEALTHCHECK --interval=30s --timeout=8s --start-period=45s --retries=3 \
  CMD ["R", "-q", "-s", "-e", "q(status = if (inherits(try(readLines('http://127.0.0.1:3838/ecm/', n = 1L, warn = FALSE), silent = TRUE), 'try-error')) 1L else 0L)"]

CMD ["/usr/local/bin/docker-entrypoint.sh"]
