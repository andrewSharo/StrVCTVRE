# This Dockerfile builds a container image for StrVCTVRE, a tool for structural variant classification.
# It uses micromamba to manage dependencies and sets up the environment for running StrVCTVRE.
# It's also been built for use with Docker and Apptainer, so it can be run on HPC clusters.
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

# The console script (`StrVCTVRE = "StrVCTVRE:main"`) imports `main` from the
# module, but the committed StrVCTVRE.py runs argparse + all logic at module
# level and defines no `main()`. Wrap the entire script in `def main():` and add
# a `__main__` guard, in-place in the image so the repo source stays untouched.
RUN awk '\
BEGIN { print "def main():" } \
{ if ($0 ~ /^[[:space:]]*$/) print $0; else print "    " $0 } \
END { print ""; print "if __name__ == \"__main__\":"; print "    main()" } \
' /software/strvctvre/StrVCTVRE.py > /tmp/StrVCTVRE.py \
  && mv /tmp/StrVCTVRE.py /software/strvctvre/StrVCTVRE.py

# Update CERTs if building on VPN
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

# Install StrVCTVRE in editable mode so that it's available as a console script (e.g. `StrVCTVRE --help` vs `python StrVCTVRE.py --help`)
RUN --mount=type=cache,target=/software/conda/pkgs \
    micromamba run -n $ENV_NAME pip --trusted-host pypi.org --trusted-host files.pythonhosted.org install -e .

WORKDIR /software/
#############################################################
#############################################################
#############################################################
# Add the micromamba activation script to the system-wide bashrc so that the environment is activated for all users on run/exec for apptainer and docker.
RUN echo "source /usr/local/bin/_activate_current_env.sh" >> /etc/bash.bashrc && \
    echo "source /usr/local/bin/_activate_current_env.sh" >> .bash_env && \
    echo "source /etc/bash.bashrc" >> .bash_env

# any calls to bash will load this environment (e.g. docker/apptainer exec bash -c "cmd")
ENV BASH_ENV=/software/.bash_env

RUN chmod +x /software/strvctvre/strvctvre_entrypoint.sh

ENV PATH=/software/strvctvre:/software:$PATH

# The `pip install -e .` above installs the `StrVCTVRE` console script, but
# setuptools flat-layout auto-discovery picks up the `data/` dir as the package
# instead of the top-level StrVCTVRE.py module, so `from StrVCTVRE import main`
# fails at runtime. Put the repo on PYTHONPATH so the module is importable.
ENV PYTHONPATH=/software/strvctvre

ENTRYPOINT ["/software/strvctvre/strvctvre_entrypoint.sh"]

CMD ["/bin/bash"]