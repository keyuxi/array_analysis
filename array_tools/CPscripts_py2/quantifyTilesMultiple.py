#!/usr/bin/env python

# Wrapper for parallelized use of AnalyseImage.m for items in a single folder
# ---------------------------------------------
#
# This script takes sequence data from an Illumina sequencer (GA or MiSeq) and registers it to image
# data collected on a GA or imaging station.
#
# Inputs:
#   CPseq files split by tile in a directory
#   Registration offset map for each tile that was imaged (.roff files; one file per tile)
#   Fluorescence Images (.tif Greyscale images; typically one image per Illumina "tile")
#
# Outputs:
#   Quantification file for each tile that was imaged (.CPfluor files; one file per tile)
#
# Curtis Layton (curtis.layton@stanford.edu) & Peter McMahon (pmcmahon@stanford.edu)
# December 2013, January 2014
# edited by Lauren Chircus (lchircus@stanford.edu) 
# July 2014, September 2014

import sys
import os
import time
import re
import argparse
import subprocess
import multiprocessing
import shutil
import uuid

import CPlibs

def collectLogs(inLog): #multiprocessing callback function to collect the output of the worker processes
    logFilename = inLog[0]
    logText = inLog[1]
    resultList[logFilename] = logText


### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description="Wrapper for parallelized use of AnalyseImage.m for items in a single folder")

parser.add_argument('-rd','--read_dir', help='directory path to the .fastq sequence read files',required=True)
parser.add_argument('-id','--image_dir', help='directory path to the tile images (default is the same as read directory)')
parser.add_argument('-ftd','--filtered_tile_dir', help='directory that holds (or will hold) the filtered sequence data (default: "read_dir/filtered_tiles")')
parser.add_argument('-rod','--roff_dir', help='directory that holds (or will hold) the registration offset files (default: "image_dir/roff")')
parser.add_argument('-n','--num_cores', help='maximum number of cores to use',required=True)
parser.add_argument('-rs','--reg_subset', help='filter names that will be used for registration (call this argument multiple times to add multiple filters to the list of filtered subsets that will be used)', action='append', required=True)
parser.add_argument('-sf','--data_scaling', help='data scaling factor.  Either "MiSeq_to_ImagingStation", "MiSeq_to_TIRFStation1", "GA_to_GA", or "none"', required=True)
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file for the run can be found', required=True)
parser.add_argument('-fd','--fluor_dir', help='directory that holds (or will hold) the CPfluor quantified fluorescence output files (default: "image_dir/CPfluor")')
parser.add_argument('-s','--stamp', help='if enabled, incorporate image stamp (text after the last underscore in the tile image filename) into the CPfluor name (default: false)', action='store_true', default=False)

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()


################ Check input files ################


#find all image files and put into a dict
if args.image_dir is not None:
    imageDir = args.image_dir
else:
    imageDir = os.path.join(args.read_dir,'images') #default image directory

if not os.path.isdir(imageDir):
    print 'Image directory "' + imageDir + '" was not found.  Please create the directory and copy tile images into it'
    sys.exit()

print 'Finding tile images in directory "' + imageDir + '"...'
tileImageDict = CPlibs.findTileFilesInDirectory(imageDir, ['.tif','.tiff'], ['.mask.tif','.mask.tiff'])

if len(tileImageDict) == 0:
    print 'No tile images files were found in the image directory "' + imageDir + '". Exiting...'
    sys.exit()

    

################ Make Filtered Tile Dictionary ################

if args.filtered_tile_dir is not None:
    filtered_tile_dir = args.filtered_tile_dir
else:
    filtered_tile_dir = os.path.join(args.read_dir,'filtered_tiles')
    
if not os.path.isdir(filtered_tile_dir):
    print 'Filtered tile directory "' + filtered_tile_dir + '" was not found.'
    sys.exit()

print 'Finding filtered tiles in "' + filtered_tile_dir + '"...'

#send out each tile to matlab to be alignment filtered (multiprocess)
filteredCPseqFilenameDict = {}   
if os.path.isdir(filtered_tile_dir): #if the merged CPseq file is not already present
    print 'Filtered tile directory "' + filtered_tile_dir + '" already exists...'
    print '   finding filtered tile .CPseq files in directory "' + filtered_tile_dir + '"...'
    filteredCPseqFilenameDict = CPlibs.findTileFilesInDirectory(filtered_tile_dir, ['.CPseq'], [])
    



################ Prepare for multicore processing ################

#assign the number of cores we will use for multicore processing based on the command-line parameter that was passed in and the number of tiles 
numCores = int(args.num_cores)
numTiles = len(filteredCPseqFilenameDict)
if (numTiles < numCores): #no need to use more cores than we actually have tiles to process
    numCores = numTiles


################ Generate the registration offset map directory ################

if args.roff_dir is not None:
    roff_dir = args.roff_dir
else:
    roff_dir = os.path.join(imageDir,'roff')
    
if not os.path.isdir(roff_dir):
    print 'Registration offset map directory "' + roff_dir + '" was not found.'
    sys.exit()

print 'Finding registration offset maps...'

roffFilenameDict = {}
maskFilenameDict = {}
if os.path.isdir(roff_dir):
    print 'Registration offset directory "' + roff_dir + '" already exists...'
    print '   finding registration offset files in directory "' + roff_dir + '"...'
    roffFilenameDict = CPlibs.findTileFilesInDirectory(roff_dir, ['.roff'], [])
    print '   finding field of view mask images in directory "' + imageDir + '"...'
    maskFilenameDict = CPlibs.findTileFilesInDirectory(imageDir, ['.mask.tif','.mask.tiff'], [])
else:
    print 'Registration offset map directory "' + roff_dir + '" was not found.'
    sys.exit()
    
    

################ Quantify cluster fluorescence from the image ################

if args.fluor_dir is not None:
    CPfluor_dir = args.fluor_dir
else:
    CPfluor_dir = os.path.join(imageDir,'CPfluor')

print 'Quantifying fluorescence of tiles...'

CPfluorFilenameDict = {}
if os.path.isdir(CPfluor_dir):
    print 'Quantified fluorescence directory "' + CPfluor_dir + '" already exists...'
    print '   finding .CPfluor files in directory "' + CPfluor_dir + '"...'
    CPfluorFilenameDict = CPlibs.findTileFilesInDirectory(CPfluor_dir, ['.CPfluor'], [])
    
    #check if any of the files found are temporary files and remove them
    for currTile, currFilename in CPfluorFilenameDict.items():
        if CPlibs.isTemporaryFilename(currFilename):
            CPfluorFilenameDict.pop(currTile,None) #remove temp file from the dict
            try:
                os.unlink(currFilename) #remove the temp file
            except OSError:
                pass
else:
    os.makedirs(CPfluor_dir) #create a new output directory if it doesn't exist

if len(CPfluorFilenameDict) < len(roffFilenameDict): #if some tiles have not been successfully quantified, attempt to quanitfy
    print '   proceeding with fluorescence quantification...'
    
    #send out each filtered tile to be quantified
    resultList = {}
    workerPool = multiprocessing.Pool(processes=numCores) #create a multiprocessing pool that uses at most the specified number of cores
    for currTileNum,currTileFilename in filteredCPseqFilenameDict.items():
        
        currCPseqFilename = filteredCPseqFilenameDict[currTileNum]
        
        if not currTileNum in tileImageDict:
            print '*** no image found for tile ' + currTileNum + '. This tile will not be quantified. ***'
            continue
        currTileImageFilename = tileImageDict[currTileNum]
        
        if not currTileNum in roffFilenameDict:
            print '*** no registration offset map found for tile ' + currTileNum + '. This tile will not be quantified. ***'
            continue
        currRoffFilename = roffFilenameDict[currTileNum]

        if args.stamp:
            currCPfluorFilename = CPlibs.getQuantifiedFilenameFromCPseqFilenameAndImageStamp(currCPseqFilename,currTileImageFilename,CPfluor_dir)
        else:
            currCPfluorFilename = CPlibs.getQuantifiedFilenameFromCPseqFilename(currCPseqFilename,CPfluor_dir)
        
        if (currTileNum in CPfluorFilenameDict) and os.path.isfile(currCPfluorFilename):
            print '*** Fluorescence quantification has already been done for tile ' + currTileNum + ' in file ' + currCPfluorFilename + '. Skipping...***'
            continue

        print 'quantifying the sequencing data in "' + currTileFilename + '"  with image "' + currTileImageFilename + '".'
        
        if currTileNum in maskFilenameDict.keys():
            currMaskFilename = maskFilenameDict[currTileNum]
            print '   using mask image "' + currMaskFilename + '".'
        else:
            currMaskFilename = ''
            print '   no mask image will be used in this quantification.'
        
        CPfluorFilenameDict[currTileNum] = currCPfluorFilename #store the list of quantified filenames

        workerPool.apply_async(CPlibs.quantifyFluorescence, args=(currCPseqFilename, currTileImageFilename, args.data_scaling, args.reg_subset, currRoffFilename, currMaskFilename, args.global_vars_path, CPlibs.tempFilename(currCPfluorFilename)), callback=collectLogs) #send out the jobs
    
    workerPool.close()
    workerPool.join()
    
    #print logs
    for currFile,currLog in resultList.items():
        print '[++++++++++++++++++++++++++++++ MATLAB ANALYSE IMAGE LOG FOR ' + currFile + ' ++++++++++++++++++++++++++++++]'
        print currLog

    #check output files
        print 'Quantified fluorescence files output to directory: "' + CPfluor_dir + '".'
    atLeastOne = False
    for currTile, currFilename in CPfluorFilenameDict.items():
        try:
            os.rename(CPlibs.tempFilename(currFilename),currFilename) #rename temp files to final filename to indicate that we are finished
        except OSError:
            pass
        
        if not os.path.isfile(currFilename): #check to make sure output files were generated
            print '      tile ' + currTile + ' file "' + currFilename + '" was NOT successfully generated.'
            CPfluorFilenameDict.pop(currTile,None) #remove it from the dict
        else:
            print '      tile ' + currTile + ' to file: "' + currFilename + '"'
            atLeastOne = True
    if atLeastOne == False:
        print 'No output files were successfully generated! Exiting...'
        sys.exit()

