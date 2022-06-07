#! /usr/bin/bash

make clean 
make 
mpirun -np 4 main.exe

