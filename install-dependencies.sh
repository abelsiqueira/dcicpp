#!/bin/bash

install_dir=$(pwd)/..
cd $install_dir
#openblas
wget http://github.com/xianyi/OpenBLAS/tarball/v0.2.9 -O openblas.tar.gz
#cholmod
wget http://www.cise.ufl.edu/research/sparse/cholmod/current/CHOLMOD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/amd/current/AMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/camd/current/CAMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/colamd/current/COLAMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/ccolamd/current/CCOLAMD.tar.gz
wget http://www.cise.ufl.edu/research/sparse/SuiteSparse_config/current/SuiteSparse_config.tar.gz
#metis
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
#others
git clone https://github.com/abelsiqueira/base_matrices.git
git clone https://github.com/abelsiqueira/nope.git

# Uncompress
for i in *.tar.gz; do tar -zxf $i; done
mv xianyi-OpenBLAS-* openblas
mv metis-4.0.3 metis

# Install
# openblas
cd openblas
make
sudo make PREFIX=/usr install
cd ..

#metis
cd metis
make
sudo cp libmetis.a /usr/lib/
cd ..

#cholmod
cd SuiteSparse_config
sed -i 's/BLAS = -lblas -lgfortran/BLAS = -lopenblas -lgfortran -lgfortranbegin -lpthread/g' SuiteSparse_config.mk
sed -i '/LAPACK = /d' SuiteSparse_config.mk
sed -i 's:METIS_PATH =.*:METIS_PATH = $install_dir:g' SuiteSparse_config.mk
sed -i 's/METIS =.*/METIS = -lmetis/g' SuiteSparse_config.mk
cd ..
for dir in CHOLMOD AMD CAMD COLAMD CCOLAMD SuiteSparse_config
do
  cd $dir
  make all
  sudo make install
  cd ..
done

#base_matrices
cd base_matrices
make all
sudo make install
cd ..

#nope
cd nope
make all
sudo make install
cd ..

#cutest - todo
