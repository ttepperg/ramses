#!/bin/bash

# Executes input make file and creates a compilation log

make -f ${1} | tee compile.log
