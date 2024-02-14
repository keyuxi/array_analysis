#!/usr/bin/env python

# Wrapper for GenerateRegistrationOffsetMap.m for multiple files in a directory
# ---------------------------------------------
#
# This script takes sequence data from an Illumina sequencer (GA or MiSeq) and registers it to image
# data collected on a GA or imaging station.
#
# Inputs:
#   CPseq files split by tile in a directory
#   Fluorescence Images (.tif Greyscale images; typically one image per Illumina "tile")
#
# Outputs:
#   Registration offset map for each cluster that was imaged (.roff files; one file per tile)
#
# Curtis Layton (curtis.layton@stanford.edu) & Peter McMahon (pmcmahon@stanford.edu)
# December 2013, January 2014
# edited by Lauren Chircus (lchircus@stanford.edu)
# July 2014, September 2014

import sys
import os
import argparse
import multiprocessing
import CPlibs

def collectLogs(inLog):
    '''
    multiprocessing callback function to collect the output of the worker processes
    '''

    logFilename = inLog[0]
    logText = inLog[1]
    resultList[logFilename] = logText

### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description="Wrapper for GenerateRegistrationOffsetMap.m for multiple files in a directory and paralelizes it")

parser.add_argument('-rd','--read_dir', help='directory path to the .fastq sequence read files')
parser.add_argument('-id','--image_dir', help='directory path to the tile images (default is the same as read directory)')
parser.add_argument('-ftd','--filtered_tile_dir', help='directory that holds (or will hold) the filtered sequence data (default: "read_dir/filtered_tiles")')
parser.add_argument('-rod','--roff_dir', help='directory that holds (or will hold) the registration offset files (default: "image_dir/roff")')
parser.add_argument('-n','--num_cores', help='maximum number of cores to use',required=True)
parser.add_argument('-rs','--reg_subset', help='filter names that will be used for registration (call this argument multiple times to add multiple filters to the list of filtered subsets that will be used)', action='append', required=True)
parser.add_argument('-sf','--data_scaling', help='data scaling factor.  Either "MiSeq_to_ImagingStation", "MiSeq_to_TIRFStation1", "GA_to_GA", or "none"', required=True)
#parser.add_argument('-sd','--flowcell_side', help='The side of the flowcell that was imaged.  Either "top" or "bottom"', required=True)
parser.add_argument('-gv','--global_vars_path', help='path to the directory in which the "GlobalVars.m" parameter file for the run can be found', required=True)

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
    print('Image directory "' + imageDir + '" was not found.  Please create the directory and copy tile images into it')
    sys.exit()

print('Finding tile images in directory "' + imageDir + '"...')
tileImageDict = CPlibs.findTileFilesInDirectory(imageDir, ['.tif','.tiff'], ['.mask.tif','.mask.tiff'])

if len(tileImageDict) == 0:
    print('No tile images files were found in the image directory "' + imageDir + '". Exiting...')
    sys.exit()


################ Make Filtered Tile Dictionary ################

if args.filtered_tile_dir is not None:
    filtered_tile_dir = args.filtered_tile_dir
else:
    filtered_tile_dir = os.path.join(args.read_dir,'filtered_tiles')

print('Finding filtered tiles in "' + filtered_tile_dir + '"...')

#send out each tile to matlab to be alignment filtered (multiprocess)
filteredCPseqFilenameDict = {}
if os.path.isdir(filtered_tile_dir): #if the merged CPseq file is not already present
    print('Filtered tile directory "' + filtered_tile_dir + '" already exists...')
    print('   finding filtered tile .CPseq files in directory "' + filtered_tile_dir + '"...')
    filteredCPseqFilenameDict = CPlibs.findTileFilesInDirectory(filtered_tile_dir, ['.CPseq'], [])




################ Prepare for multicore processing ################

#assign the number of cores we will use for multicore processing based on the command-line parameter that was passed in and the number of tiles
numCores = int(args.num_cores)
numTiles = len(filteredCPseqFilenameDict)
if (numTiles < numCores): #no need to use more cores than we actually have tiles to process
    numCores = numTiles


################ Generate the registration offset maps ################

if args.roff_dir is not None:
    roff_dir = args.roff_dir
else:
    roff_dir = os.path.join(imageDir,'roff')

print('Generating registration offset maps from filtered tile-specific sequences and tile image files...')

roffFilenameDict = {}
maskFilenameDict = {}
if os.path.isdir(roff_dir):
    print('Registration offset directory "' + roff_dir + '" already exists...')
    print('   finding registration offset files in directory "' + roff_dir + '"...')
    roffFilenameDict = CPlibs.findTileFilesInDirectory(roff_dir, ['.roff'], [])
    print('   finding field of view mask images in directory "' + imageDir + '"...')
    maskFilenameDict = CPlibs.findTileFilesInDirectory(imageDir, ['.mask.tif','.mask.tiff'], [])
else:
    os.makedirs(roff_dir) #create a new output directory if it doesn't exist

if len(roffFilenameDict) < len(filteredCPseqFilenameDict): #if some or all registration offset maps have not yet been generated, attempt to generate
    print('   proceeding with registration offset map generation...')

    #generate registration offset maps and mask images
    resultList = {}
    workerPool = multiprocessing.Pool(processes=numCores) #create a multiprocessing pool that uses at most the specified number of cores
    for currTileNum,currTileFilename in list(filteredCPseqFilenameDict.items()):

        currCPseqFilename = filteredCPseqFilenameDict[currTileNum]

        if not currTileNum in tileImageDict:
            print('*** no image found for tile ' + currTileNum + '. A registration map will not be generated. ***')
            continue
        currTileImageFilename = tileImageDict[currTileNum]

        currRoffFilename = CPlibs.getRoffFilenameFromImageFilename(currTileImageFilename,roff_dir)
        currMaskFilename = CPlibs.getMaskFilenameFromImageFilename(currTileImageFilename) #mask files (if generated) are stored in the image directory

        if (currTileNum in roffFilenameDict) and os.path.isfile(currRoffFilename):
            print('*** a registration  offset map has already been generated for tile ' + currTileNum + ' in file ' + currRoffFilename + '. Skipping...***')
            continue

        print('generating the registration offset map for "' + currTileFilename + '"  and image "' + currTileImageFilename + '".')

        roffFilenameDict[currTileNum] = currRoffFilename #store the list of registration offset map filenames
        maskFilenameDict[currTileNum] = currMaskFilename #store the list of field of view mask images

        workerPool.apply_async(CPlibs.getRegistrationOffset, args=(currCPseqFilename, currTileImageFilename, args.data_scaling, args.reg_subset, currRoffFilename, currMaskFilename, args.global_vars_path), callback=collectLogs) #send out the job

    workerPool.close()
    workerPool.join()

    #print logs
    for currFile,currLog in list(resultList.items()):
        print('[++++++++++++++++++++++++++++++ MATLAB GENERATE REGISTRATION OFFSET MAP LOG FOR ' + currFile + ' ++++++++++++++++++++++++++++++]')
        print(currLog)

    #check output files
    print('Registration offset map files output to directory: "' + roff_dir + '".')
    atLeastOne = False
    for currTile,currFilename in list(roffFilenameDict.items()):
        if not os.path.isfile(currFilename): #check to make sure output files were generated
            print('      tile ' + currTile + ' file "' + currFilename + '" was NOT successfully generated.')
            roffFilenameDict.pop(currTile,None) #remove it from the dict
        else:
            print('      tile ' + currTile + ' to file: "' + currFilename + '"')
            atLeastOne = True
    if atLeastOne == False:
        print('No output files were successfully generated! Exiting...')
        sys.exit()

    print('If generated, mask images are in the image directory "' + imageDir + '".')

