@echo off
setlocal
git config --global core.longpaths true
whoami

mkdir build
mkdir build/release
mkdir build/debug


cd cmake
conan install . --profile=conan_profile_win64_debug --remote=lab --output-folder=../build/debug --build=missing --update
conan install . --profile=conan_profile_win64_debug --remote=lab --output-folder=../build/release --build=missing --update
cd ..