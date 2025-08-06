#!/bin/bash

# Creates a list of my standard patched files, excluding dice, dice_restart, oagertz, and restart.
# This list is mostly relevant to the patch oagertz. See oagertz/README.

prefix='mypatches_list'

# delete previous list
rm ${prefix}.txt

ls ../{amr,hydro,pm,poisson}/*f90 | awk -F'[/]' '{print $(NF)}' > ${prefix}.tmp

# sort alphabetically and remove duplicates (shouldn't be any)
sort -u ${prefix}.tmp > ${prefix}.txt

echo -e "\n\tThese are the standard patched files:\n"

cat ${prefix}.txt

# delete temporary list
rm ${prefix}.tmp

echo -e "\n\tThis list is contained in ${prefix}.txt"

echo -e "\tDone.\n"
