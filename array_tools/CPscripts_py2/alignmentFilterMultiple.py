#!/usr/bin/env python

# Wrapper for alignmentFilter() function: runs MATLAB instances that perform alignment filters on every CPseq file in a given directory

import sys
import os
import argparse
import multiprocessing
import CPlibs



def collectLogs(inLog): #multiprocessing callback function to collect the output of the worker processes
    logFilename = inLog[0]
    logText = inLog[1]
    resultList[logFilename] = logText


### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description="Wrapper for alignmentFilter() function: runs MATLAB instances that perform alignment filters on every CPseq file in a given directory")

parser.add_argument('-rd','--read_dir', help='directory in which CPseq data files to filter are located',required=True)
parser.add_argument('-f','--filter', help='filename of the .CPfilter file',required=True)
parser.add_argument('-od','--output_dir', help='directory to write filtered CPseq data files to',required=True)
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file can be found', required=True)
parser.add_argument('-n','--num_cores', help='maximum number of threads to use simultaneously', required=True)

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()

numCores = int(args.num_cores)

#check input directory
if not os.path.isdir(args.read_dir):
    print 'input directory "' + args.read_dir + '" was not found.  Exiting...'
    sys.exit()

#check filter file
if not os.path.isfile(args.filter):
    print '.CPfilter file "' + args.filter + '" was not found.  Exiting...'
    sys.exit()

#check if output directory exists
if not os.path.isdir(args.output_dir):
    # make new directory
    os.makedirs(args.output_dir) 

print 'Performing alignment-based filtering based on sequences and rules from file "' + args.filter + '"...'

#send out each file to matlab to be alignment filtered (multiprocess)
    
resultList = {} #dict to hold results
workerPool = multiprocessing.Pool(processes=numCores) #create a multiprocessing pool that uses at most the specified number of cores
for f in os.listdir(args.read_dir): # use listdir() instead of walk(), so that it doesn't dig into subdirectories
    if os.path.isfile(os.path.join(args.read_dir,f)):
        fileExt = os.path.splitext(f)[-1] # Split the extension from the path
        if fileExt.lower() == '.cpseq': 
            currFullFilenameInput = os.path.join(args.read_dir, f)
            currFullFilenameOutput = CPlibs.getFilteredFilenameFromTileFilename(f,args.output_dir)
            workerPool.apply_async(CPlibs.alignmentFilter, args=(currFullFilenameInput, args.filter, currFullFilenameOutput, args.global_vars_path), callback=collectLogs) #send out the job

workerPool.close()
workerPool.join()

#print logs
for currFile,currLog in resultList.items():
    print '[++++++++++++++++++++++++++++++ MATLAB ALIGNMENT LOG FOR ' + currFile + ' ++++++++++++++++++++++++++++++]'
    print currLog
