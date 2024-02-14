#!/usr/bin/env python

# Wrapper for alignmentFilter() function: runs a MATLAB instance that performs an alignment filter on a single CPseq file

import sys
import os
import argparse
import CPlibs

### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description="Wrapper for alignmentFilter() function: runs a MATLAB instance that performs an alignment filter on a single CPseq file")

parser.add_argument('-i','--input', help='filename of CPseq data to filter',required=True)
parser.add_argument('-f','--filter', help='filename of the .CPfilter file',required=True)
parser.add_argument('-o','--output', help='filename to write filtered CPseq data to',required=True)
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file can be found', required=True)

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()

#check input CPseq file
if not os.path.isfile(args.input):
    print 'input .CPseq file "' + args.input + '" was not found.  Exiting...'
    sys.exit()

#check filter file
if not os.path.isfile(args.filter):
    print '.CPfilter file "' + args.filter + '" was not found.  Exiting...'
    sys.exit()

logString = CPlibs.alignmentFilter(args.input, args.filter, args.output, args.global_vars_path)
print(logString)