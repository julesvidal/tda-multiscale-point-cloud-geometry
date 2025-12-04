#! /bin/bash

mkdir install
cd src/
mkdir build
cd build/

cmake -DCMAKE_INSTALL_PREFIX=../../install ..
make -j install

cd ..

