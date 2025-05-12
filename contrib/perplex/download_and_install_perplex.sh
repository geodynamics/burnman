#!/bin/bash

if [ -d perplex-installer ]; then
    cd perplex-installer
    git pull --recurse-submodules origin main
    cd ..
else
    git clone --recurse-submodules https://github.com/bobmyhill/perplex-installer.git
    cd perplex-installer
    ./install_perplex.sh
    cd ..
fi
