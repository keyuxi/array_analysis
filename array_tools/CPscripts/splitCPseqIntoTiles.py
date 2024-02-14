#!/usr/bin/env python

#written by Curtis Layton Nov 2013
#splits a CPseq file (or other file where the first token of every line is a cluster id) into an idividual file for every tile
#the cluster id (first field of the tab-delimited file) should be of the format "<machine id>:<run index>:<flowcell id>:<lane #>:<tile id>:<x coord>:<y coord>"

import os
import sys
import argparse
import re

import CPlibs

def getNumLinesInFile(filename):
    n_lines = 0
    with open(filename) as f:
        for line in f: n_lines += 1
    return n_lines

def update_progress(currSequenceIdx, numSequences):
    updateInterval = 50
    if(currSequenceIdx % updateInterval == 0) | (currSequenceIdx+1 == numSequences):
        barLength = 50 # Modify this to change the length of the progress bar
        progress = float(currSequenceIdx+1)/float(numSequences)
        status = ""
        if currSequenceIdx+1 == numSequences:
            status = "Done...\r\n"
        block = int(round(barLength*progress))
        text = "\rPercent: [{0}] {1:3.2f}% ({2} sequences) {3}".format( "#"*block + "-"*(barLength-block), progress*100, currSequenceIdx+1, status)
        sys.stdout.write(text)
        sys.stdout.flush()

def checkValidOutputFilename(outputFilename):
    with open(outputFilename,"w") as f: #attempt to open for writing.  will raise error if fails
        pass
    try:
        os.unlink(outputFilename) #remove the file that was temporarily created
    except OSError:
        pass

def getOutputFilename(inFullPath,currTile,currSide,outputPath):
    #split the input filename into parts
    (inPath,inFilename) = os.path.split(inFullPath)
    (inBasename,ext) = os.path.splitext(inFilename)
    
    if currSide == 1:
        currSideString = 'Top'
    elif currSide == 2:
        currSideString = 'Bottom'
    
    outputFilename = os.path.join(outputPath,inBasename + '_tile' + currTile + '_' + currSideString + ext) #append tile designation
    return outputFilename

def outputToTileFile(currInputLine, currTile, currSide, outputFilenames, inputFilename, outputPath):
    tileId = currTile + str(currSide) #note if the swath is to be considered, this will have to change
    if (not outputFilenames) or (tileId not in outputFilenames):
        #if the filename for this tile does not already exist, create it
        outputFilenames[tileId] = getOutputFilename(inputFilename,currTile,currSide,outputPath)
        with open(CPlibs.tempFilename(outputFilenames[tileId]),"w") as outputFileHandle: #create file with temporary filename indicating the file is incomplete
            outputFileHandle.write(currInputLine)
    else:
        with open(CPlibs.tempFilename(outputFilenames[tileId]),"a") as outputFileHandle: #append to temporary file
            outputFileHandle.write(currInputLine)
    return outputFilenames




### MAIN ###
def main(rawArgs, **keywordParams):
    
    #set up command line argument parser
    parser = argparse.ArgumentParser(description="Splits a .CPseq file into individual .CPseq files for each tile. Output filenames are derived automatically from the input filename.")
    
    parser.add_argument('inputFilename', help='input .CPseq file to be split into tile-specific .CPseq files')
    parser.add_argument('-o','--output_dir', help='output directory. Default is the same directory as the input file.')
    parser.add_argument('-s','--side', help='specify if tiles will be output from the top surface (default), bottom surface, or both.', choices=['top','bottom','both'], default='top')
    outputExtension = '.CPseq'
    
    #parse command line arguments
    args = parser.parse_args(rawArgs[1:])

    #get the number of lines in the input file
    numLinesInInputFile = getNumLinesInFile(args.inputFilename)
    
    #get output directory
    if args.output_dir is None:
        (outputPath, inputFilenameStripped) = os.path.split(args.inputFilename) #use the same path as the input file
    else:
        outputPath = args.output_dir
    
    #create a dictionary to hold output filenames
    outputFilenames = {}
    
    #open the input file and begin reading
    print('Processing ' + str(numLinesInInputFile) + ' entries from input file: ' + args.inputFilename)
    with open(args.inputFilename) as inputFile:
        progressCount = 0;
        for currInputLine in inputFile:
            
            currID = re.split('[\t ]',currInputLine)[0] #get the cluster ID (first token of the line)

            (flowcellSN,currSide,currSwath,currTile) = CPlibs.parseClusterID(currID)

            #NOTE: the swath is ignored in this implementation.
            
            currTile = '{:03}'.format(int(currTile)) #format currTile in 3-digit string in anticipation of porting the pipeline to sytems with more tiles
            
            if args.side.lower() == 'both':
                outputFilenames = outputToTileFile(currInputLine,currTile,currSide,outputFilenames,args.inputFilename,outputPath)
            elif(args.side.lower() == 'top' and currSide == 1):
                outputFilenames = outputToTileFile(currInputLine,currTile,currSide,outputFilenames,args.inputFilename,outputPath)
            elif(args.side.lower() == 'bottom' and currSide == 2):
                outputFilenames = outputToTileFile(currInputLine,currTile,currSide,outputFilenames,args.inputFilename,outputPath)
                
            update_progress(progressCount, numLinesInInputFile)
            progressCount = progressCount + 1
    
    #rename all files to final filename indicating that the files are complete
    for currFilename in list(outputFilenames.values()):
        os.rename(CPlibs.tempFilename(currFilename),currFilename)
    
    if ('from_command_line' in keywordParams):
        #function was called from the command line
        return 0 #return exit value
    else:
        #function was imported and called from a library
        return outputFilenames #return filename dictionary of all output files generated

if __name__=='__main__':
    sys.exit(main(sys.argv, from_command_line=True))


