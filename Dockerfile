FROM mambaorg/micromamba:ubuntu24.04 AS micromamba

FROM ubuntu:24.04

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    wget \
    curl \
    ca-certificates \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -f /var/cache/apt/archives/*.deb

#############################################################
######################### General ###########################
#############################################################
RUN mkdir -p /software; chmod a+rwx -R /software; mkdir /software/strvctvre

WORKDIR /software/strvctvre

COPY . .

# Update CERTs if on VPN
# RUN wget -P /tmp/ <cert_url>
# ARG CURL_CA_BUNDLE=/tmp/<cert_file>

# Install micromamba first
# Adding micromamba to existing image:  https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html#adding-micromamba-to-an-existing-docker-image
# https://github.com/mamba-org/micromamba-docker/blob/main/debian.Dockerfile

# if your image defaults to a non-root user, then you may want to make the
# next 3 ARG commands match the values in your image. You can get the values
# by running: docker run --rm -it my/image id -a
ARG MAMBA_USER=root
ARG MAMBA_USER_ID=0
ARG MAMBA_USER_GID=0
ENV MAMBA_USER=$MAMBA_USER
ENV MAMBA_ROOT_PREFIX="/software/conda"
ENV MAMBA_EXE="/bin/micromamba"
## conda env name (must match `name:` in environment_py3.yml)
ENV ENV_NAME=StrVCTVRE_py_3
#ENV MAMBA_SKIP_ACTIVATE=1; We want to activate, but just not the base

COPY --from=micromamba "$MAMBA_EXE" "$MAMBA_EXE"
COPY --from=micromamba /usr/local/bin/_activate_current_env.sh /usr/local/bin/_activate_current_env.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_shell.sh /usr/local/bin/_dockerfile_shell.sh
COPY --from=micromamba /usr/local/bin/_entrypoint.sh /usr/local/bin/_entrypoint.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN /usr/local/bin/_dockerfile_initialize_user_accounts.sh && \
    /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN --mount=type=cache,target=/software/conda/pkgs \
    micromamba create -y -f /software/strvctvre/environment_py3.yml && \
    micromamba clean --all --yes

#############################################################
#############################################################
#############################################################

ENV PATH=/software/strvctvre:/software:$PATH

RUN chmod +x /software/strvctvre/strvctvre_entrypoint.sh

ENTRYPOINT ["/software/strvctvre/strvctvre_entrypoint.sh"]