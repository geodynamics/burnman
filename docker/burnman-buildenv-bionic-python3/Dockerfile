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
    texlive \
    texlive-latex-extra
