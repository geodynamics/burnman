FROM ubuntu:18.04

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
    numdiff \
    python3 \
    python3-matplotlib \
    python3-numpy \
    python3-scipy \
    python3-sympy \
    python3-pip \
    libgmp3-dev \
    texlive \
    texlive-latex-extra

# The following are not available as an Ubuntu package. Note that pycddlib
# requires libgmp3-dev.

RUN pip3 install pycddlib
