#!/bin/bash

#clean


#run conan
mkdir build
cd build
mkdir release
mkdir debug
cd ../cmake
/home/gitlab-runner/anaconda3/envs/conan/bin/conan install . --profile=conan_profile_linux_debug --remote=lab --output-folder=../build/debug --build=missing --update
/home/gitlab-runner/anaconda3/envs/conan/bin/conan install . --profile=conan_profile_linux_debug --remote=lab --output-folder=../build/release --build=missing --update