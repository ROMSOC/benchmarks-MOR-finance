# gnu-octave/docker

# Please follow docker best practices
# https://docs.docker.com/engine/userguide/eng-image/dockerfile_best-practices/

## Synchronized with https://github.com/jupyter/docker-stacks/tree/master/base-notebook

FROM  gnuoctave/octave:6.4.0
LABEL maintainer="Kai T. Ohlhus <k.ohlhus@gmail.com>"

ENV LAST_UPDATED=2022-02-01

ARG NB_USER="jovyan"
ARG NB_UID="1000"
ARG NB_GID="100"

# Fix DL4006
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install JupyterLab

USER root

RUN chmod 777 /tmp        && \
    apt-get update --yes  && \
    DEBIAN_FRONTEND="noninteractive" \
    apt-get install --yes --no-install-recommends \
      imagemagick            \
      tini                && \
    # Install required python packages
    #   - sympy, see https://savannah.gnu.org/bugs/?58491
    pip3 install --upgrade --no-cache-dir \
      jupyterlab                    \
      octave_kernel                 \
      jupytext                      \
      jupyter-book                  \
      ghp-import                    \
      numpy                         \
      sympy==1.5.1                  \
      matplotlib                 && \
    apt-get --yes clean          && \
    apt-get --yes autoremove     && \
    rm -Rf /var/lib/apt/lists/*

# Configure environment
ENV SHELL=/bin/bash \
    NB_USER="${NB_USER}" \
    NB_UID=${NB_UID} \
    NB_GID=${NB_GID} \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    HOME="/home/${NB_USER}"

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
# Switch back to jovyan to avoid accidental container runs as root
USER ${NB_UID}

WORKDIR "${HOME}"
