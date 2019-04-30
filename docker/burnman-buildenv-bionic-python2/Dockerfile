FROM ubuntu:18.04

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
    numdiff \
    python \
    python-matplotlib \
    python-numpy \
    python-scipy \
    python-sympy \
    texlive \
    texlive-latex-extra
