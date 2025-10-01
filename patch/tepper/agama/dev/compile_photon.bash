#!/bin/bash

# Link against the dynamic library and compile
# Arguments are source name and exec name, in that order!
# The 'l:agama.so' syntax means agama.so literally and not libagama.so
#
# !!!!!!!!!! IMPORTANT !!!!!!!!!!
# Make sure to add AGAMA_PATH to LD_LIBRARY_PATH
# BASH:
# export LD_LIBRARY_PATH=${AGAMA_PATH}:${LD_LIBRARY_PATH}
# CSH:
# setenv LD_LIBRARY_PATH ${AGAMA_PATH}:${LD_LIBRARY_PATH}
# Note: the rpath linking flag (as used below) makes this unnecessary

FC='gfortran'
AGAMA_PATH='/suphys/tepper/.local/lib/python3.9/site-packages/agama'
LIB_FLAGS='-Wl,-rpath,'
LIBAGAMA='-L'${AGAMA_PATH}' '${LIB_FLAGS}${AGAMA_PATH}'  -l:agama.so'
# The following may or may not be needed:
#LIBOMP='-fopenmp'

if test $# -ne 2 ; then
    echo -e '\n\tNeed to provide f90 filename and exec filename, in that order!\n'
  exit
fi

echo -e '\n\tRemoving previous exec '${2}
rm ${2}

command=${FC}' -o '${2}' '${1}' '${LIBAGAMA}' '${LIBOMP}

echo -e '\n\tCompiling with: \n\t'${command}'\n'
echo ${command} | sh
