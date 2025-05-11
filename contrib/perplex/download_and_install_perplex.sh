#!/bin/bash

if [ -d perplex-installer ]; then
    rm -rf perplex-installer
fi
git clone --recursive git@github.com:bobmyhill/perplex-installer.git
cd perplex-installer
./install_perplex.sh
cd ..