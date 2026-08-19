#!/bin/bash

# Cleans after executing a make file with make_log.bash

make -f ${1} clean

# in case compilation failed, this might be present
rm -f write_makefile.f90
rm -f write_patch.f90

# remove customo compile log
rmlog='rm -f compile.log'
echo ${rmlog}
echo ${rmlog} | sh
