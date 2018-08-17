#!/bin/bash

# Change this to match your compiler
CXX="clang++"
OPT="-O3 -g"
FEATURES="-std=c++14 -fpic"
WARNS="-Wall -Wno-sign-compare"
PYTHON="python3.6"


${PYTHON} GetCoreCount.py
${PYTHON} CreateGaussLegendreTable.py
${CXX} -lpthread -shared -o libgyro.so ${FEATURES} ${OPT} ${WARNS} Simulation.cpp
