#!/bin/bash

# Automatically retrieve the source files required for this patch.
# Note that these files are the ORIGINAL source files from
# Oscar Agertz' Ramses patch and they NEED TO BE MODIFIED to
# work with my patches as outlined in README_NEW

# The retrieval is done in two ways:
# 1 - Manual
# 2 - Semi-automatic

# Run: $> ./retrieve_src_orig.sh
# Make sure it's executable ($> chmod +x retrieve_src.bash)

# delete previous list
rm src_orig_list.txt

# 1 - Manual

# Needed by feedback.f90
echo '/Users/tepper/codes/ramses_agertz/ramses/pm/metal_yields.f90' >> src_orig_list.tmp
# Not needed, but convenient:
echo '/Users/tepper/codes/ramses_agertz/ramses/amr/movie.f90' >> src_orig_list.tmp
# Needed if movie.f90 is added:
echo '/Users/tepper/codes/ramses_agertz/ramses/pm/particle_snapshot.f90' >> src_orig_list.tmp
# Needed by particle_snapshot.f90 (included if movie.f90 is added)
echo '/Users/tepper/codes/ramses_agertz/ramses/pm/tracer_utils.f90' >> src_orig_list.tmp


# 2 - Semi-automatic

# Explanation: Relevant files have the additional, non-standard variables: mpb (particle birth mass), 'zp' (particle metallicity), and additional metal variable names (identifiable mey_key); these are flagged alternatively by 'EDGE' or 'ERIC', i.e. the simulation suite and developer who introduced them to Ramses, or the non-standard variable 'nmetals'. Moreover, there are some tweaks related to tracer variables; these are flagged by 'tracer_ivar' or 'tracer_var'. Finally, a custom-made log file relevant to star-formation and supernovae events are flagged by 'SFunit' and 'SNunit', respectively. Therefore, a viable approach to retrieve the relevant source files is to look for those containing any of these strings.

# Set common look-up path
DIRS='/Users/tepper/codes/ramses_agertz/ramses/{amr,hydro,io,pm,poisson}'

# Info to stdout
echo -e '\n\tLooking for files in '$DIRS'\n'

# grep + regex: looks for substring 'mpb' avoding letter before or after
# 1 awk: filters out the first column of the output
# 2 awk: extracts substring terminated by ':', effectively the filename with full path
substr='mpb'
echo "grep -rE '(^|[^[:alpha:]])${substr}([^[:alpha:]]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

substr='ERIC'
echo "grep -rE '(^|[^[:alpha:]])${substr}([^[:alpha:]]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

substr='EDGE'
echo "grep -rE '(^|[^[:alpha:]])${substr}([^[:alpha:]]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

substr='met_keys'
echo "grep -rE '(^|[^[:alpha:]])${substr}([^[:alpha:]]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

substr='nmetals'
echo "grep -rE '(^|[^[:alpha:]])${substr}([^[:alpha:]]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

# grep + regex: looks for substring 'zp' avoding '_' after
substr='zp'
echo "grep -rE '(^|[^[:alpha:]])${substr}([^[:alpha:]_]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

# grep + regex: looks for substring 'tracer_var' or 'tracer_ivar'
substr='S[NF]unit'
echo "grep -rE 'tracer_(i)?var' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

# grep + regex: looks for substring 'SFunit' or 'SNunit' avoding letter before or after
substr='S[NF]unit'
echo "grep -rE '(^|[^[:alpha:]])S[NF]unit([^[:alpha:]]|$)' ${DIRS}" | sh | awk '{print $1}' | awk -F'[:]' '{print $(NF-1)}' >> src_orig_list.tmp

# sort file and remove duplicates
sort -u src_orig_list.tmp > src_orig_list.txt

echo -e "\n\tFound the following files:\n"
cat src_orig_list.txt | awk '{print "\t"$1}'

# RSYNC files to destination dir
dest='src_orig'
echo -e '\n\tRSYNCing files to destination dir '${dest}
command='rsync -av --files-from=src_orig_list.txt --no-relative / '${dest}'/.'
echo -e '\n\t'${command}
echo -e '\n\t'${command} | sh

# delete temporary list
rm src_orig_list.tmp

echo -e "\n\tDone.\n"
