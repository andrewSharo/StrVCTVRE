#!/bin/bash

# StrVCTVRE container entrypoint.
# Activates the conda environment, then warns if the required phyloP
# conservation file is missing (it is ~9.2GB and is NOT baked into the image).

IGreen='\033[0;92m'       # Green
On_Red='\033[41m'         # Red
Color_Off='\033[0m'       # Text Reset

# Path StrVCTVRE.py checks by default; override at runtime with -p.
PHYLOP_FILE="${PHYLOP_FILE:-data/hg38.phyloP100way.bw}"

if [ ! -s "${PHYLOP_FILE}" ]; then
  echo -e "${IGreen}The phyloP conservation file was not found at:\n  ${PHYLOP_FILE}${Color_Off}"
  echo -e "${IGreen}StrVCTVRE requires hg38.phyloP100way.bw (~9.2GB). Download it from UCSC:${Color_Off}"
  echo -e "${IGreen}  https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw${Color_Off}"
  echo -e "${On_Red}\tMount it into the container (e.g. docker: -v /host/data; apptainer: --bind /host/data)${Color_Off}"
  echo -e "${On_Red}\tor pass its path to StrVCTVRE.py with -p.${Color_Off}"
fi

# Run whatever command was passed to the container (e.g. python StrVCTVRE.py ...),
# or drop into an interactive shell if none was given.
if [ "$#" -gt 0 ]; then
  source /usr/local/bin/_activate_current_env.sh
  exec "$@"
else
  exec /bin/bash
fi
