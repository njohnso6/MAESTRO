#!/bin/bash

# install MAESTRO/python
$PYTHON setup.py install

# install giggle
cd refpkg/giggle
make
cp bin/giggle $PREFIX/bin/
cd ../..

# install sinto, the pypi version is not useful, let's
# do it through git
git clone https://github.com/timoast/sinto
cd sinto
git checkout 16b5336 # for version 0.7.2
$PYTHON setup.py install
cd ../

# there are two dependencies in R DESCRIPTION
# that can't be found in conda-forge or bioconda
# channel. They are grid and Gmisc

# install MAESTRO/R
$R -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true");devtools::install(".", upgrade="never")'


