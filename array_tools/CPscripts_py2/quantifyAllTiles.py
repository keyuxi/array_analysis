#!/usr/bin/env python

# Main Pipeline for Fluorescence Quantification
# ---------------------------------------------
#
# This high-level script takes sequence data from an Illumina sequencer (GA or MiSeq) and matches sequences
# from those data with quantified fluorescence values of the clusters on the same flowcell.
#
# Inputs:
#   Sequence data (.fastq files)
#   Filter (.CPfilter file; used to define subsets of the sequence data, typically based on sequence)
#   Fluorescence Images (.tif Greyscale images; typically one image per Illumina "tile")
#
# Outputs:
#   Quantified fluorescence value for each cluster that was imaged (.CPfluor files; one file per tile)
#
# There are a number of intermediate files that are produced in the process of going from the initial input
# to the final output. This script will "pick up where it left off" if a failure occurs inbetween the start
# and end of the whole operation.
# 
# Similarly, it is possible to use this script to perform quantification of several sets of images taken of
# the same flowcell (e.g. at different temperatures, or concentrations) without needing to re-run steps that
# are independent of the images (e.g. conversion of .fastq files to .CPseq files; filtering of the CPseq files;
# potentially the generation of registration offset maps [.roff files]).
#
#
# Curtis Layton (curtis.layton@stanford.edu) & Peter McMahon (pmcmahon@stanford.edu)
# December 2013, January 2014

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


#top-level master script for running the image analysis pipeline

### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description="master script for the CP image analysis pipeline")

parser.add_argument('-rd','--read_dir', help='directory path to the .fastq sequence read files',required=True)
parser.add_argument('-td','--tile_dir', help='directory that holds (or will hold) CPseq data from each tile (default: "read_dir/tiles")')
parser.add_argument('-id','--image_dir', help='directory path to the tile images (default is the same as read directory)')
parser.add_argument('-f','--filter', help='filename of the .CPfilter file',required=True)
parser.add_argument('-ftd','--filtered_tile_dir', help='directory that holds (or will hold) the filtered sequence data (default: "read_dir/filtered_tiles")')
parser.add_argument('-rod','--roff_dir', help='directory that holds (or will hold) the registration offset files (default: "image_dir/roff")')
parser.add_argument('-fd','--fluor_dir', help='directory that holds (or will hold) the CPfluor quantified fluorescence output files (default: "image_dir/CPfluor")')
parser.add_argument('-n','--num_cores', help='maximum number of cores to use',required=True)
parser.add_argument('-rs','--reg_subset', help='filter names that will be used for registration (call this argument multiple times to add multiple filters to the list of filtered subsets that will be used)', action='append', required=True)
parser.add_argument('-sf','--data_scaling', help='data scaling factor.  Either "MiSeq_to_ImagingStation", "MiSeq_to_TIRFStation1", "GA_to_GA", or "none"', required=True)
parser.add_argument('-sd','--flowcell_side', help='The side of the flowcell that was imaged.  Either "top" or "bottom"', required=True)
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file for the run can be found', required=True)
parser.add_argument('-s','--stamp', help='if enabled, incorporate image stamp (text after the last underscore in the tile image filename) into the CPfluor name (default: false)', action='store_true', default=False)

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()


################ Check input files ################

#check filter file
if not os.path.isfile(args.filter):
    print '.CPfilter file "' + args.filter + '" was not found.  Exiting...'
    sys.exit()

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



################ Merge fastqs into CPseq format ################

#get output filename from the flowcell ID in the input sequence (e.g. fastq) files
inputFilenames = CPlibs.findInputSequenceFilesInDirectory(args.read_dir,[[],[],[],[]],silentMode=True,requiredToFind=True);
firstFilename = CPlibs.getFirstFilenameInInputFilenames(inputFilenames)
(masterFlowcellID,inputFilePath) = CPlibs.getFlowcellSNfromFile(firstFilename,'fastq') #get Flowcell ID

mergedCPseqFilename = os.path.join(args.read_dir,masterFlowcellID + '_ALL.CPseq') #name merged CPseq file according to flowcell ID

print 'Merging fastq files in directory "' + args.read_dir + '" into file: "' + mergedCPseqFilename + '"...'

if os.path.isfile(mergedCPseqFilename): #if the merged CPseq file is not already present
    print 'File "' + mergedCPseqFilename + '" already exists. Skipping...'
else:
    if os.path.isfile(CPlibs.tempFilename(mergedCPseqFilename)): #delete temp file if already exists
        os.unlink(CPlibs.tempFilename(mergedCPseqFilename))
    CPlibs.mergeFastqs(args.read_dir,CPlibs.tempFilename(mergedCPseqFilename))
    try:
        os.rename(CPlibs.tempFilename(mergedCPseqFilename),mergedCPseqFilename) #rename to indicate that it finished
    except OSError:
        print 'merged .CPseq file "' + mergedCPseqFilename + '" was not successfully generated.  Exiting...'
        sys.exit()



################ Split into individual tile CPseqs ################

if args.tile_dir is not None:
    tile_dir = args.tile_dir
else:
    tile_dir = os.path.join(args.read_dir,'tiles')

print 'Splitting CPseq file "' + mergedCPseqFilename + '" into indivicual tile files...'

tempFilesFound = False
tileCPseqFilenameDict = {}
if os.path.isdir(tile_dir): #if the merged CPseq file is already present
    print 'Tile directory "' + tile_dir + '" already exists...'
    print '   finding tile .CPseq files in directory "' + tile_dir + '"...'
    tileCPseqFilenameDict = CPlibs.findTileFilesInDirectory(tile_dir, ['.CPseq'], [])

    #check if any of the files found are temporary files and remove them
    for currTile, currFilename in tileCPseqFilenameDict.items():
        if CPlibs.isTemporaryFilename(currFilename):
            tempFilesFound = True
            tileCPseqFilenameDict.pop(currTile,None) #remove temp file from the dict
            try:
                os.unlink(currFilename) #remove the temp file
            except OSError:
                pass
else:
    os.makedirs(tile_dir) #create a new output directory if it doesn't exist

if len(tileCPseqFilenameDict) == 0 or tempFilesFound:
    print '   proceeding to split into tiles...'
    tileCPseqFilenameDict = CPlibs.splitCPseq(args.flowcell_side, mergedCPseqFilename, tile_dir) #split into tiles in the temp directory
    print '   Individual tile files output to directory: "' + tile_dir + '".'
else:
    print '   Individual tile files already generated.  Skipping...  Final tile files:'
    
#check output files
atLeastOne = False
for currTile, currFilename in tileCPseqFilenameDict.items():
    if not os.path.isfile(currFilename): #check to make sure output files were generated
        print '      tile ' + currTile + ' file "' + currFilename + '" was NOT successfully generated.'
        tileCPseqFilenameDict.pop(currTile,None) #remove it from the dict
    else:
        print '      tile ' + currTile + ' to file: "' + currFilename + '"'
        atLeastOne = True
if atLeastOne == False:
    print 'No tile .CPseq files were successfully generated! Exiting...'
    sys.exit()

################ Prepare for multicore processing ################

#assign the number of cores we will use for multicore processing based on the command-line parameter that was passed in and the number of tiles 
numCores = int(args.num_cores)
numTiles = len(tileCPseqFilenameDict)
if (numTiles < numCores): #no need to use more cores than we actually have tiles to process
    numCores = numTiles

################ Do alignment-based sequence filtering ################

if args.filtered_tile_dir is not None:
    filtered_tile_dir = args.filtered_tile_dir
else:
    filtered_tile_dir = os.path.join(args.read_dir,'filtered_tiles')

print 'Performing alignment-based filtering based on sequences and rules from file "' + args.filter + '"...'

#send out each tile to matlab to be alignment filtered (multiprocess)
filteredCPseqFilenameDict = {}   
if os.path.isdir(filtered_tile_dir): #if the merged CPseq file is not already present
    print 'Filtered tile directory "' + filtered_tile_dir + '" already exists...'
    print '   finding filtered tile .CPseq files in directory "' + filtered_tile_dir + '"...'
    filteredCPseqFilenameDict = CPlibs.findTileFilesInDirectory(filtered_tile_dir, ['_filtered.CPseq'], [])
    
    #check if any of the files found are temporary files and remove them
    for currTile, currFilename in filteredCPseqFilenameDict.items():
        if CPlibs.isTemporaryFilename(currFilename):
            filteredCPseqFilenameDict.pop(currTile,None) #remove temp file from the dict
            try:
                os.unlink(currFilename) #remove the temp file
            except OSError:
                pass
else:
    os.makedirs(filtered_tile_dir) #create a new output directory if it doesn't exist

if len(filteredCPseqFilenameDict) < len(tileCPseqFilenameDict): #if some or all filtered .CPseq files have not been generated, attempt to generate them
    print '   proceeding with alignment filtering...'
    
    resultList = {} #dict to hold results
    workerPool = multiprocessing.Pool(processes=numCores) #create a multiprocessing pool that uses at most the specified number of cores
    for currTileNum,currTileFilename in tileCPseqFilenameDict.items():
        
        currFilteredFilename = CPlibs.getFilteredFilenameFromTileFilename(currTileFilename,filtered_tile_dir)
        
        if (currTileNum in filteredCPseqFilenameDict) and os.path.isfile(currFilteredFilename):
            print '*** tile ' + currTileNum + ' has already been alignment-filtered to file ' + currFilteredFilename + '. Skipping...***'
            continue
        
        filteredCPseqFilenameDict[currTileNum] = currFilteredFilename #store the list of filtered filenames
        
        workerPool.apply_async(CPlibs.alignmentFilter, args=(currTileFilename, args.filter, CPlibs.tempFilename(currFilteredFilename), args.global_vars_path), callback=collectLogs) #send out the job

    workerPool.close()
    workerPool.join()

    #print logs
    for currFile,currLog in resultList.items():
        print '[++++++++++++++++++++++++++++++ MATLAB ALIGNMENT LOG FOR ' + currFile + ' ++++++++++++++++++++++++++++++]'
        print currLog

    #check output files
    print 'Alignment filtered tile files output to  directory: "' + filtered_tile_dir + '".'
    atLeastOne = False
    for currTile, currFilename in filteredCPseqFilenameDict.items():
        try:
            os.rename(CPlibs.tempFilename(currFilename),currFilename) #rename temp files to final filename to indicate that we are finished
        except OSError:
            pass
        
        if not os.path.isfile(currFilename): #check to make sure output files were generated
            print '      tile ' + currTile + ' file "' + currFilename + '" was NOT successfully generated.'
            filteredCPseqFilenameDict.pop(currTile,None) #remove it from the dict
        else:
            print '      tile ' + currTile + ' to file: "' + currFilename + '"'
            atLeastOne = True
    if atLeastOne == False:
        print 'No output files were successfully generated! Exiting...'
        sys.exit()

################ Generate the registration offset maps ################

if args.roff_dir is not None:
    roff_dir = args.roff_dir
else:
    roff_dir = os.path.join(imageDir,'roff')

print 'Generating registration offset maps from filtered tile-specific sequences and tile image files...'

roffFilenameDict = {}
maskFilenameDict = {}
if os.path.isdir(roff_dir):
    print 'Registration offset directory "' + roff_dir + '" already exists...'
    print '   finding registration offset files in directory "' + roff_dir + '"...'
    roffFilenameDict = CPlibs.findTileFilesInDirectory(roff_dir, ['.roff'], [])
    print '   finding field of view mask images in directory "' + imageDir + '"...'
    maskFilenameDict = CPlibs.findTileFilesInDirectory(imageDir, ['.mask.tif','.mask.tiff'], [])
else:
    os.makedirs(roff_dir) #create a new output directory if it doesn't exist

if len(roffFilenameDict) < len(filteredCPseqFilenameDict): #if some or all registration offset maps have not yet been generated, attempt to generate
    print '   proceeding with registration offset map generation...'
    
    #generate registration offset maps and mask images
    resultList = {}
    workerPool = multiprocessing.Pool(processes=numCores) #create a multiprocessing pool that uses at most the specified number of cores
    for currTileNum,currTileFilename in filteredCPseqFilenameDict.items():
        
        currCPseqFilename = filteredCPseqFilenameDict[currTileNum]
        
        if not currTileNum in tileImageDict:
            print '*** no image found for tile ' + currTileNum + '. A registration map will not be generated. ***'
            continue
        currTileImageFilename = tileImageDict[currTileNum]

        currRoffFilename = CPlibs.getRoffFilenameFromImageFilename(currTileImageFilename,roff_dir)
        currMaskFilename = CPlibs.getMaskFilenameFromImageFilename(currTileImageFilename) #mask files (if generated) are stored in the image directory
        
        if (currTileNum in roffFilenameDict) and os.path.isfile(currRoffFilename):
            print '*** a registration  offset map has already been generated for tile ' + currTileNum + ' in file ' + currRoffFilename + '. Skipping...***'
            continue

        print 'generating the registration offset map for "' + currTileFilename + '"  and image "' + currTileImageFilename + '".'
         
        roffFilenameDict[currTileNum] = currRoffFilename #store the list of registration offset map filenames
        maskFilenameDict[currTileNum] = currMaskFilename #store the list of field of view mask images
               
        workerPool.apply_async(CPlibs.getRegistrationOffset, args=(currCPseqFilename, currTileImageFilename, args.data_scaling, args.reg_subset, currRoffFilename, currMaskFilename, args.global_vars_path), callback=collectLogs) #send out the job
            
    workerPool.close()
    workerPool.join()
    
    #print logs
    for currFile,currLog in resultList.items():
        print '[++++++++++++++++++++++++++++++ MATLAB GENERATE REGISTRATION OFFSET MAP LOG FOR ' + currFile + ' ++++++++++++++++++++++++++++++]'
        print currLog

    #check output files
    print 'Registration offset map files output to directory: "' + roff_dir + '".'
    atLeastOne = False
    for currTile,currFilename in roffFilenameDict.items():
        if not os.path.isfile(currFilename): #check to make sure output files were generated
            print '      tile ' + currTile + ' file "' + currFilename + '" was NOT successfully generated.'
            roffFilenameDict.pop(currTile,None) #remove it from the dict
        else:
            print '      tile ' + currTile + ' to file: "' + currFilename + '"'
            atLeastOne = True
    if atLeastOne == False:
        print 'No output files were successfully generated! Exiting...'
        sys.exit()
        
    print 'If generated, mask images are in the image directory "' + imageDir + '".'

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
