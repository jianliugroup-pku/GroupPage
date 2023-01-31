#!/bin/bash 
gfortran -O2 EOMs_Liu.f90 -llapack -o eom_L
gfortran -O2 EOMs_Huo1.f90 -llapack -o eom_H1
gfortran -O2 EOMs_Huo2.f90 -llapack -o eom_H2
gfortran -O2 EOMs_Liu_CI.f90 -llapack -o eom_CI

# run code examples
./eom_L weak-SW
./eom_L strong-SW
./eom_H1 weak-SW-Huo1
./eom_H1 strong-SW-Huo1
./eom_H2 weak-SW-Huo2
./eom_H2 strong-SW-Huo2
./eom_CI Berry-SW 

# plot result
python3 plot.py
