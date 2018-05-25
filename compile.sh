#!/bin/bash

# Change this to match your compiler
CXX="clang++"
OPT="-O3 -g"
FEATURES="-std=c++14 -fpic"
WARNS="-Wall -Wno-sign-compare"

${CXX} -lpthread -shared -o libgyro.so ${FEATURES} ${OPT} ${WARNS} Simulation.cpp
