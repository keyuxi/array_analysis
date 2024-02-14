#!/usr/bin/env python

# Pipeline for Fluorescence Quantification with previously generated CPseq files
# ---------------------------------------------
#
# This script will run the tail end of the 'quantifyAllTiles.py' script. 
# It requires that filtered CPseq files have already been created.
#
# Inputs:
#   Sequence data (.CPseq files)
#   Fluorescence Images (.tif Greyscale images; still one image per Illumina "tile")
#
# Outputs:
#   Quantified fluorescence value for each cluster that was imaged (.CPfluor files; one file per tile)
#
# 
#
# Curtis Layton (curtis.layton@stanford.edu) & Peter McMahon (pmcmahon@stanford.edu)
# and now Ben Ober-Reynolds (boberrey@stanford.edu)
# December 2013, January 2014, August 2016

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

parser.add_argument('-id','--image_dir', help='directory path to the tile images.', required=True)
parser.add_argument('-ftd','--filtered_tile_dir', help='directory that holds the filtered sequence data.', required=True)
parser.add_argument('-rod','--roff_dir', help='directory that holds (or will hold) the registration offset files (default: "image_dir/roff")')
parser.add_argument('-fd','--fluor_dir', help='directory that holds (or will hold) the CPfluor quantified fluorescence output files (default: "image_dir/CPfluor")')
parser.add_argument('-n','--num_cores', help='maximum number of cores to use',required=True)
parser.add_argument('-rs','--reg_subset', help='filter names that will be used for registration (call this argument multiple times to add multiple filters to the list of filtered subsets that will be used)', action='append', required=True)
parser.add_argument('-sf','--data_scaling', help='data scaling factor.  Either "MiSeq_to_ImagingStation", "MiSeq_to_TIRFStation1", "GA_to_GA", or "none"', required=True)
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file for the run can be found', required=True)

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()

#parse command line arguments
args = parser.parse_args()


################ Check input files ################


imageDir = args.image_dir
if not os.path.isdir(imageDir):
    print 'Image directory "' + imageDir + '" was not found.  Please create the directory and copy tile images into it'
    sys.exit()

#find all image files and put into a dict

print 'Finding tile images in directory "' + imageDir + '"...'
tileImageDict = CPlibs.findTileFilesInDirectory(imageDir, ['.tif','.tiff'], ['.mask.tif','.mask.tiff'])

if len(tileImageDict) == 0:
    print 'No tile images files were found in the image directory "' + imageDir + '". Exiting...'
    sys.exit()


# find all filtered CPseq files and put into a dict

filtered_tile_dir = args.filtered_tile_dir
if not os.path.isdir(filtered_tile_dir):
    print 'filtered_tile_dir directory "' + filtered_tile_dir + '" was not found.  Please create the directory and copy filtered CPseq files into it'
    sys.exit()

print 'Finding CPseq files in directory "' + filtered_tile_dir + '"...'
filteredCPseqFilenameDict = CPlibs.findTileFilesInDirectory(filtered_tile_dir, ['.CPseq'], [])

if len(tileImageDict) == 0:
    print 'No CPseq files were found in the filtered_tile_dir directory "' + filtered_tile_dir + '". Exiting...'
    sys.exit()




################ Generate the registration offset maps ################

# prepare for multiprocessing:
numCores = int(args.num_cores)
numTiles = len(filteredCPseqFilenameDict)
if (numTiles < numCores): #no need to use more cores than we actually have tiles to process
    numCores = numTiles


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
        currCPfluorFilename = CPlibs.getQuantifiedFilenameFromCPseqFilename(currTileImageFilename, CPfluor_dir)
        
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
