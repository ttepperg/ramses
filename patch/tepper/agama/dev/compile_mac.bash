#!/bin/bash

# Link against the dynamic library and compile
# Arguments are source name and exec name , in that order!
#
# !!!!!! IMPORTANT !!!!!!
# Make sure to add AGAMA_PATH to DYLD_LIBRARY_PATH:
# BASH:
# export DYLD_LIBRARY_PATH=${AGAMA_PATH}:${DYLD_LIBRARY_PATH}
# CSH:
# setenv DYLD_LIBRARY_PATH ${AGAMA_PATH}:${DYLD_LIBRARY_PATH}

if test $# -ne 2 ; then
  echo -e '\n\tNeed to provide f90 filename and exec filename, in that order!\n'
  exit
fi

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo -e "\n\tMUST source this script, i.e. run as:
            $> source compile_mac.bash <f90> <exe>\n"
    exit
fi

AGAMA_PATH='/Users/tepper/Library/Python/3.9/lib/python/site-packages/agama'

# Check if DIR is already in LD_LIBRARY_PATH
if [[ ":$DYLD_LIBRARY_PATH:" == *":$AGAMA_PATH:"* ]]; then
    echo -e "\n\t$AGAMA_PATH is already in DYLD_LIBRARY_PATH"
else
    echo -e "\n\tAdding $AGAMA_PATH to DYLD_LIBRARY_PATH"
    export DYLD_LIBRARY_PATH="$AGAMA_PATH:${DYLD_LIBRARY_PATH:-}"
fi

FC='gfortran-mp-14'
LIBAGAMA='-L'${AGAMA_PATH}' -Wl,-rpath,'${AGAMA_PATH}' '${AGAMA_PATH}'/agama.so'
# The following may or may not be needed:
# LIBOMP='-L/opt/local/lib/libomp -lomp'
LIBOMP=

echo -e '\n\tRemoving previous exec '${2}
rm ${2}

command=${FC}' -o '${2}' '${1}' '${LIBAGAMA}' '${LIBOMP}

echo -e '\n\tCompiling with: \n\t'${command}'\n'
echo ${command} | sh
