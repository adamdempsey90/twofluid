#!/bin/bash
echo "Compiling "$1".c"
gcc-4.9 -c -Wall src/$1.c -o bin/$1.o

