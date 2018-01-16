#!/bin/bash

mkdir bin

g++ -std=c++11 -Iexternal/cinolib/include -Iexternal/eigen3 main.cpp -obin/hyper_surface_smoothing