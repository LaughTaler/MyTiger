#!/bin/bash
scl enable devtoolset-7 bash

cd build
cd release
cmake -DCMAKE_INSTALL_PREFIX=~/tiger/release ../.. -Wno-dev -DCMAKE_BUILD_TYPE=Release
make -j 
make install
#
cd ../debug
cmake -DCMAKE_INSTALL_PREFIX=~/tiger/debug ../.. -Wno-dev -DCMAKE_BUILD_TYPE=Debug
make -j 
make install

