#!/bin/bash

ALL=$1 	#~ a:all, s:single
FORT=$2 	#~ gfortran, gfortran-7, etc.
DEBUG=$3	#~ yes/no
GPROF=$4	#~ yes/no

if [ $ALL == "a" ]
then
	make clean -f makefile
fi

make "FORT=$2" "DEBUG=$3" "GPROF=$4" -f makefile
