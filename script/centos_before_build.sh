#!/bin/bash
#use gcc-7
scl enable devtoolset-7 bash

#clean

git reset --hard HEAD
rm src/blmesh -rf
rm src/meshgen -rf
rm src/wyfdt -rf
rm src/DT_remesh -rf
rm src/extrude -rf
rm src/meshmetric -rf
rm src/trans2_geo -rf
rm src/Mesh_Repair -rf

#run conan
mkdir build
cd build
mkdir release
mkdir debug
cd ../cmake

/home/fuwuqi1/anaconda3/envs/conan/bin/conan install . --profile=conan_profile_linux_debug --remote=lab --output-folder=../build/debug --build=missing
/home/fuwuqi1/anaconda3/envs/conan/bin/conan install . --profile=conan_profile_linux_debug --remote=lab --output-folder=../build/release --build=missing
