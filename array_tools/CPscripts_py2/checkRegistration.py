#!/usr/bin/env python

# Wrapper for checkTileRegistration.m 
# ---------------------------------------------
#
# for quickly checking registration of a single image to the .CPseq sequence data for a tile
# Curtis Layton Jan 2015

import sys
import os
import argparse
import CPlibs
from cStringIO import StringIO

#create a string that defines a matlab cell array from a python array of strings
def assembleMatlabCellArray(pythonArrayOfStrings):
    print pythonArrayOfStrings
    output = "{";
    if pythonArrayOfStrings:
        for i, currStr in enumerate(pythonArrayOfStrings):
            output = output + "'" + currStr + "'"
            if i < (len(pythonArrayOfStrings)-1): #except for the last item..
                output = output + "," #add a comma
    output = output + "}";
    return output

### MAIN ###

#unique string that is matched to checkTileRegistration.m that is used to identify that the program ran normally and produced the expected output 
checkRegUUID = '9d7bd2bb-5f26-4893-89ac-641e30ec8ab9';

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description="Wrapper for checkTileRegistration.m for quickly checking registration of a single image to the .CPseq sequence data for a tile")

parser.add_argument('-s','--seq_filename', help='filename of the sequence data for the tile (.CPseq format)',required=True)
parser.add_argument('-i','--image_filename', help='filename of the tile image to be registered', required=True)
parser.add_argument('-sf','--data_scaling', help='data scaling factor', choices=['MiSeq_to_ImagingStation', 'MiSeq_to_TIRFStation1', 'GA_to_GA', 'none'], required=True)
parser.add_argument('-f','--filter_subset', nargs='+', help='filter subsets to be used for image registration. These filters should be present in the .CPseq file. ')
parser.add_argument('-o','--output_image_path', help='path for the output images (if not specified images will not be generated)')
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file for the run can be found', required=True)
parser.add_argument('-r','--remote', help='formats the output to be machine readable for remote calls', action='store_true')

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()


################ Check input params ################

if not os.path.isfile(args.seq_filename):
    print 'ERROR: sequence tile .CPseq file "' + args.seq_filename + '" was not found.  Aborting'
    sys.exit()

if not os.path.isfile(args.image_filename):
    print 'ERROR: Image file "' + args.image_filename + '" was not found.  Aborting'
    sys.exit()

if args.output_image_path is not None:
    if not os.path.isdir(args.output_image_path):
        print 'ERROR: Output image path "' + args.output_image_path + '" does not exist. Aborting'
        sys.exit()

globalVarsFilename = os.path.join(args.global_vars_path, 'GlobalVars.m')
if not os.path.isfile(globalVarsFilename):
    print 'ERROR: GlobalVars.m file not found in directory "' + args.global_vars_pat + '".  Aborting'
    sys.exit()
    
################ Call checkTileRegistration.m in matlab ################

#redirect stdio into a string to supress the output of the job spawning, etc.
sys.stdout = redirectedStdout = StringIO()

#assemble the filter subset identifiers in matlab cell array format
filterSubsetString = assembleMatlabCellArray(args.filter_subset)

#format output image path
if args.output_image_path is None:
    outputImagePath = [];
else:
    outputImagePath = "'" + args.output_image_path + "'"
    
#assemble the matlab function call
if args.remote:
    matlabFunctionCallString = "checkTileRegistration('{0}','{1}','{2}', {3}, {4},true);".format(args.seq_filename, args.image_filename, args.data_scaling, filterSubsetString, outputImagePath)
else:
    matlabFunctionCallString = "checkTileRegistration('{0}','{1}','{2}', {3}, {4});".format(args.seq_filename, args.image_filename, args.data_scaling, filterSubsetString, outputImagePath)

#submit the job to be run by matlab 
matlabOutput = CPlibs.spawnMatlabJob(matlabFunctionCallString, args.global_vars_path);

#direct stdout back to normal for final output
sys.stdout = sys.__stdout__

#format output
if args.remote:
    #split the output based on a unique ID output by checkTileRegistration.m
    splitMatlabOutput = matlabOutput.rsplit(checkRegUUID,1)
    if (len(splitMatlabOutput) > 1):
        #if the output ID was recognized and used to parse the output
        print checkRegUUID + splitMatlabOutput[1]
    else:
        #if not, something went wrong--just barf out the output for debugging
        print matlabOutput
else:
    print matlabOutput

