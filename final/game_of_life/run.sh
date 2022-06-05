#! /usr/bin/bash

make clean 
make 
mpirun -np 5 main.exe

