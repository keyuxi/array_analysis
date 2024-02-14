#!/usr/bin/env python

# written by Curtis Layton Nov 2013 to merge multiple fastq files that came FROM THE SAME RUN (e.g. a paired-end run)
# this is intended to consolidate reads that pertain to the same cluster into one file, i.e. read1, read2, index read1, index read2 with the same cluster id.
# Files are written into a new format,".CPseq", which is a tab-delimited flat file of the format:
#
# <cluster ID> <filter IDs> <read 1> <phred quality read 1> <read 2> <phred quality read 2> <index read 1> <phred quality index read 1> <index read 2> <phred quality index read 2>
#
# where the format of <cluster ID> is:
# <machine id>:<run index>:<flowcell id>:<lane #>:<tile #>:<x coord>:<y coord>
#
# The program colates sequences by cluster id and does not require that they be in the same order
# All clusters are outputted, including those that are present in some files but missing in others

import os
import sys
import time
import re
from Bio import SeqIO
import argparse
import shelve
import uuid
import itertools

import CPlibs

filetype = 'fastq'
numLinesPerRecord = 4 #there are 4 lines for every record in a fastq file
outputExtension = '.CPseq'
memoryThreshold = 1000000 #max number of records that we can load into memory before we use the "shelf" database for collation
tempShelfFilename = 'cluster_shelf_' + str(uuid.uuid4()) + str(time.time()) + '.dat'

def phredStr(phredArray):
    if all(x == 0 for x in phredArray):
        return ''
    else:
        phredQuality = [chr(q+33) for q in phredArray]
        return ''.join(phredQuality)

def getNumLinesInFile(filename):
    n_lines = 0
    with open(filename) as f:
        for line in f: n_lines += 1
    return n_lines

class ClusterData:
    def __init__(self):
        self.filterID = ''
        self.read1 =  ''
        self.qRead1 = []
        self.read2 = ''
        self.qRead2 = []
        self.index1 = ''
        self.qIndex1 = []
        self.index2 = ''
        self.qIndex2 = []

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

def getMaxRecordsPerRead(numRecordsInFile):
    maxRecordsPerRead = 0
    for currRead in numRecordsInFile:
        currRecordsPerRead = 0
        for currNumRecords in currRead:
            currRecordsPerRead = currRecordsPerRead + currNumRecords
        maxRecordsPerRead = max(maxRecordsPerRead, currRecordsPerRead)
    return maxRecordsPerRead

#get a filter ID from the standard Illumina filename convention
def getFilterIDfromIlluminaFilename(filename):
    (path,filename) = os.path.split(filename) #strip off the path
    matches = re.finditer('_S[0-9]*_L[0-9]*_[IR][0-9]*_[0-9]*',filename) #search for the illumina filename pattern
    indices = [m.start(0) for m in matches]
    if(indices == []):
        filterID = os.path.splitext(filename)[0]; #if no matches were found just use the whole filename (stripping off the file extension)
    else:
        filterID = filename[0:indices[-1]]; #else use the filename prefix

    if filterID.lower() == 'undetermined':
        return ''
    else:
        return filterID

def filterIDisInString(filterID, filterIDstring):
    return (filterID in filterIDstring.split(':'))

def tempFilename(fullPath): #appends a prefix to the filename to indicate the the file is incomplete
    (path,filename) = os.path.split(fullPath) #split the file into parts
    return os.path.join(path,'__' + filename)

def findMatchingFilename(inputFilename,filenameList):
    inputFilterID = getFilterIDfromIlluminaFilename(inputFilename)
    
    if os.path.isfile(inputFilename): #check if the file exists
        inputHandle = open(inputFilename,'r') #open the file for reading
        inputHandle.seek(0) #ensure we are at the beginning of the file        
        inputRecord = SeqIO.parse(inputFilename, filetype).next() #get the first record in the file
    
    for i,currFilename in enumerate(filenameList):
        
        #if the filename is in Illumina Format, get the filterID from the filename and compare to the filenames in the list to match
        currFilterID = getFilterIDfromIlluminaFilename(currFilename)
        if currFilterID == inputFilterID: #if we found a match
            return i #return the index of the matching filename
        
        #if the filename is not matched by name convention, try to open the files and compare the clusterID of the first line to match
        if os.path.isfile(currFilename): #check if the file exists
            currHandle = open(currFilename,'r') #open the file for reading
            currHandle.seek(0) #ensure we are at the beginning of the file
            currRecord = SeqIO.parse(currHandle, filetype).next()
            if currRecord.id == inputRecord.id: #compare the clusterID of the first entry in the files
                return i #return the index of the matching filename
        
    return -1 #return -1 as a flag that we did not find a matching filename.

### MAIN ###
def main(rawArgs, **keywordParams):

    #set up command line argument parser
    parser = argparse.ArgumentParser(description="Consolidates individual fastq files from a paired-end run (optionally with index reads, as well) into one file (.CPseq format)")

    parser.add_argument('-r1','--read1', help='read 1 fastq filename', action='append')
    parser.add_argument('-r2','--read2', help='read 2 fastq filename', action='append')
    parser.add_argument('-i1','--index1', help='index read 1 fastq filename', action='append')
    parser.add_argument('-i2','--index2', help='index read 2 fastq filename', action='append')
    parser.add_argument('-rd','--read_dir',help='assigns all files in the directory matching the Illumina read1 and read2 filename convention to be read1 and read2 files, respectively')
    parser.add_argument('-o','--output', help='output filename (.CPseq)')
    parser.add_argument('-m','--memory', help='force data to be read into memory (which is faster), even if above the threshold for using the disk-based shelf database',action='store_true')

    #parse command line arguments
    args = parser.parse_args(rawArgs[1:])

    if not len(sys.argv) > 1:
        parser.print_help()

    inputFilenames = [[],[],[],[]]
    numRecordsInFile = [[],[],[],[]]

    #store the input filneames in a 2D list
    if args.read1 is not None:
        inputFilenames[0] = args.read1
    if args.read2 is not None:
        inputFilenames[1] = args.read2
    if args.index1 is not None:
        inputFilenames[2] = args.index1
    if args.index2 is not None:
        inputFilenames[3] = args.index2

    #read in all files following the illumina filename convention from a specified directory
    if args.read_dir is not None:
        inputFilenames = CPlibs.findInputSequenceFilesInDirectory(args.read_dir,inputFilenames)

    if inputFilenames == [[],[],[],[]]: #if no input files were identified
        raise IOError("at least one input file must be passed in or present in the read directory")

    #get the flowcell ID
    firstFilename = CPlibs.getFirstFilenameInInputFilenames(inputFilenames)
    (masterFlowcellID,inputFilePath) = CPlibs.getFlowcellSNfromFile(firstFilename,filetype)

    #get the number of records in each file, get read lengths, and check flowcell IDs
    for i in range(len(inputFilenames[0])):
        numRecordsInFile[0].append(getNumLinesInFile(inputFilenames[0][i])/numLinesPerRecord)
        if(numRecordsInFile[0][i] > 0): #if the file is not empty
            currHandle = open(inputFilenames[0][i],'r') #open the file for reading
            currHandle.seek(0) #ensure we are at the beginning of the file   
            currRecord = SeqIO.parse(currHandle, filetype).next() #read the first record in the file
            #readLengths[0] = max(readLengths[0],len(currRecord.seq))
            if CPlibs.getFlowcellSNfromClusterID(currRecord.id).lower() != masterFlowcellID.lower():
                raise IOError('All input files to be merged together must be from the same run. File ' + inputFilenames[0][i] + ' does not match.')
    for i in range(len(inputFilenames[1])):
        numRecordsInFile[1].append(getNumLinesInFile(inputFilenames[1][i])/numLinesPerRecord)
        if(numRecordsInFile[1][i] > 0): #if the file is not empty
            currHandle = open(inputFilenames[1][i],'r') #open the file for reading
            currHandle.seek(0) #ensure we are at the beginning of the file   
            currRecord = SeqIO.parse(currHandle, filetype).next() #read the first record in the file
            #readLengths[1] = max(readLengths[1],len(currRecord.seq))
            if CPlibs.getFlowcellSNfromClusterID(currRecord.id).lower() != masterFlowcellID.lower():
                raise IOError('All input files to be merged together must be from the same run. File ' + inputFilenames[1][i] + ' does not match.')
    for i in range(len(inputFilenames[2])):
        numRecordsInFile[2].append(getNumLinesInFile(inputFilenames[2][i])/numLinesPerRecord)
        if(numRecordsInFile[2][i] > 0): #if the file is not empty
            currHandle = open(inputFilenames[2][i],'r') #open the file for reading
            currHandle.seek(0) #ensure we are at the beginning of the file   
            currRecord = SeqIO.parse(currHandle, filetype).next() #read the first record in the file
            #readLengths[2] = max(readLengths[2],len(currRecord.seq))
            if CPlibs.getFlowcellSNfromClusterID(currRecord.id).lower() != masterFlowcellID.lower():
                raise IOError('All input files to be merged together must be from the same run. File ' + inputFilenames[2][i] + ' does not match.')
    for i in range(len(inputFilenames[3])):
        numRecordsInFile[3].append(getNumLinesInFile(inputFilenames[3][i])/numLinesPerRecord)
        if(numRecordsInFile[3][i] > 0): #if the file is not empty
            currHandle = open(inputFilenames[3][i],'r') #open the file for reading
            currHandle.seek(0) #ensure we are at the beginning of the file   
            currRecord = SeqIO.parse(currHandle, filetype).next() #read the first record in the file
            #readLengths[3] = max(readLengths[3],len(currRecord.seq))
            if CPlibs.getFlowcellSNfromClusterID(currRecord.id).lower() != masterFlowcellID.lower():
                raise IOError('All input files to be merged together must be from the same run. File ' + inputFilenames[3][i] + ' does not match.')

    #check output filename
    if args.output is not None:
        currExt = os.path.splitext(args.output)[-1]
        if(currExt.lower() != outputExtension.lower()):
            raise IOError('the output file must have a \'' + outputExtension + '\' extension')
        else:
            outputFilename = args.output
    else:
        outputFilename = os.path.join(inputFilePath,masterFlowcellID + '_ALL.CPseq') #default name according to flowcell ID
    checkValidOutputFilename(outputFilename) #check that output filename is valid
    print 'Merged sequences will be saved to : ' + outputFilename



    ### attempt line-by-line colation.  This only works if all read1 and read2 files are paired and have the same clusterIDs in the exact same order in both files
    print 'attempting line-by-line collation.  This is fast and will work if all read1 and read2 files contain all the same records in the same order'
    lineByLineFailed = False

    #check that the number of reads matches between files
    #if (len(inputFilenames[0]) != len(inputFilenames[1])):
    #    print '    Number of read1 and read2 files does not match.  Quitting line-by-line collation...'
    #    lineByLineFailed = True
    #
    #if not lineByLineFailed:
    #    index1ReadsPresent = False
    #    if (len(inputFilenames[2]) != 0):
    #        if (len(inputFilenames[2]) == len(inputFilenames[0])): #check to make sure the number of index1 reads matches read1 and read2
    #            index1ReadsPresent = True
    #        else:
    #            print '    Number of index1 reads do not match read1 and read2.  Quitting line-by-line collation...'
    #            lineByLineFailed = True
    #
    #if not lineByLineFailed:
    #    index2ReadsPresent = False
    #    if (len(inputFilenames[3]) != 0):
    #        if (len(inputFilenames[3]) == len(inputFilenames[0])): #check to make sure the number of index2 reads matches read1 and read2
    #            index2ReadsPresent = True
    #        else:
    #            print '    Number of index2 reads do not match read1 and read2.  Quitting line-by-line collation...'
    #            lineByLineFailed = True


    #read in sequence files and consolidate records
    if not lineByLineFailed:
        with open(tempFilename(outputFilename),"w") as outputFileHandle: #open the output file for writing with a temporary filename to indicate not finished
            for inputFileIndex in range(len(inputFilenames[0])): #loop through all read1 files

                currRead1Filename = inputFilenames[0][inputFileIndex] #get current read1 filename
                currRead1NumRecords = numRecordsInFile[0][inputFileIndex] #get num records in the current read1 filename
                currRead1FilterID = getFilterIDfromIlluminaFilename(currRead1Filename) #get a filterID from the index-read filtered Illumina read1 filename

                dummyFilename = currRead1Filename;

                #search for the matching filename in the read2 records
                read2FileIndex = findMatchingFilename(currRead1Filename,inputFilenames[1])
                if read2FileIndex != -1: #if the file was found...
                    currRead2Filename = inputFilenames[1][read2FileIndex]
                    currRead2NumRecords = numRecordsInFile[1][read2FileIndex]
                    if currRead1NumRecords != currRead2NumRecords:
                        lineByLineFailed = True #the number of records must match for line-by-line colation
                else:
                    currRead2Filename = dummyFilename

                #search for the matching filename in the index1 records
                index1FileIndex = findMatchingFilename(currRead1Filename,inputFilenames[2])
                if index1FileIndex != -1: #if the file was found...
                    currIndex1Filename = inputFilenames[2][index1FileIndex]
                    currIndex1NumRecords = numRecordsInFile[2][index1FileIndex]
                    if currRead1NumRecords != currIndex1NumRecords:
                        lineByLineFailed = True #the number of records must match for line-by-line colation
                else:
                    currIndex1Filename = dummyFilename

                #search for the matching filename in the index2 records
                index2FileIndex = findMatchingFilename(currRead1Filename,inputFilenames[3])
                if index2FileIndex != -1: #if the file was found...
                    currIndex2Filename = inputFilenames[3][index2FileIndex]
                    currIndex2NumRecords = numRecordsInFile[3][index2FileIndex]
                    if currRead1NumRecords != currIndex2NumRecords:
                        lineByLineFailed = True #the number of records must match for line-by-line colation
                else:
                    currIndex2Filename = dummyFilename

                #if (read2FileIndex == -1) & (index1FileIndex == -1) & (index2FileIndex == -1):
                #    print '    No matching files found for read1 file "' + currRead1Filename + '".  Line-by-line collation failed.'
                #    lineByLineFailed = True

                if not lineByLineFailed:
                    print '    Colating read1 file "' + currRead1Filename + '" with:'
                    if read2FileIndex != -1:
                        print '        read2 file "' + currRead2Filename + '"'
                    if index1FileIndex != -1:
                        print '        index1 file "' + currIndex1Filename + '"'
                    if index2FileIndex != -1:
                        print '        index2 file "' + currIndex2Filename + '"'
                    print '        line-by-line into file "' + outputFilename + '"...'

                    print 'processing ' + str() + ' sequences...'
                    i = 0

                    # iterate through all files in lockstep
                    # any files that were not found are just given the read1 filename as a "dummy" file which is then just iterated through redundantly but not read from more than once
                    for (currRead1Record,currRead2Record,currIndex1Record,currIndex2Record) in itertools.izip(SeqIO.parse(currRead1Filename, filetype),SeqIO.parse(currRead2Filename, filetype),SeqIO.parse(currIndex1Filename, filetype),SeqIO.parse(currIndex2Filename, filetype)):

                        if currRead1Record.id == currRead2Record.id:
                            cl = ClusterData() #create a new cluster

                            currID = currRead1Record.id
                            cl.filterID = currRead1FilterID
                            cl.read1 = currRead1Record.seq
                            cl.qRead1 = currRead1Record.letter_annotations['phred_quality']
                            if read2FileIndex != -1:
                                cl.read2 = currRead2Record.seq
                                cl.qRead2 = currRead2Record.letter_annotations['phred_quality']
                            if index1FileIndex != -1:
                                cl.index1 = currIndex1Record.seq
                                cl.qIndex1 = currIndex1Record.letter_annotations['phred_quality']
                            if index2FileIndex != -1:
                                cl.index2 = currIndex2Record.seq
                                cl.qIndex2 = currIndex2Record.letter_annotations['phred_quality']

                            #write to file
                            outputFileHandle.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(currID, cl.filterID, cl.read1, phredStr(cl.qRead1), cl.read2, phredStr(cl.qRead2), cl.index1,phredStr(cl.qIndex1), cl.index2,phredStr(cl.qIndex2)))
                        else:
                            lineByLineFailed = True
                            print
                            print '    Sequences in file are not present in the same order.  Line-by-line collation failed.'
                            break

                        update_progress(i,currRead1NumRecords)
                        i = i + 1

                if lineByLineFailed:
                    break

        if lineByLineFailed:
            #delete temp output file
            try:
                os.unlink(CPlibs.tempFilename(outputFilename))
            except Exception:
                pass
        else:
            #line-by-line was successful, rename the temp output filename to the final filename and exit
            os.rename(CPlibs.tempFilename(outputFilename),outputFilename)
            print '    Line-by-line collation successful. Output saved to file " ' + outputFilename + '"'



    ### colate by hashing
    if lineByLineFailed:
        #if line-by-line colation failed, do the colation by hashing

        print 'colating records by hashing'

        #determine whether we can read records into memory, or if there are too many and we need to use the shelf
        maxRecordsPerRead = getMaxRecordsPerRead(numRecordsInFile)
        useDatabase = (maxRecordsPerRead > memoryThreshold) | args.memory

        if(useDatabase):
            print 'collating data on disk...'
            allClusterData = shelve.open(tempShelfFilename,'n') #init database file to hold all the data
        else:
            print 'collating data in memory...'
            allClusterData = {} #init dict to hold all data in memory

        #read in sequence files and consolidate records
        for inputReadIndex in range(len(inputFilenames)):
            for inputFileIndex in range(len(inputFilenames[inputReadIndex])):
                print 'loading file: ' + inputFilenames[inputReadIndex][inputFileIndex]
                print 'processing ' + str(numRecordsInFile[inputReadIndex][inputFileIndex]) + ' sequences...'
                i = 0

                #get a filterID from the index-read filtered Illumina filename
                currFilterID = getFilterIDfromIlluminaFilename(inputFilenames[inputReadIndex][inputFileIndex])

                for currRecord in SeqIO.parse(inputFilenames[inputReadIndex][inputFileIndex], filetype):
                    if currRecord.id in allClusterData: #if found, add it to the existing record
                        currCluster = allClusterData[currRecord.id] #get a local reference to the matching ClusterData object
                    else: #if not found, create a new cluster
                        currCluster = ClusterData() #new cluster

                    #Illumina names files according to index-read filtering. Use the filename to add a filterID
                    if(currCluster.filterID == ''):
                        currCluster.filterID = currFilterID
                    elif(filterIDisInString(currFilterID, currCluster.filterID)): #do not duplicate if the filterID is already present
                        pass
                    else:
                        currCluster.filterID = currCluster.filterID + ':' + currFilterID

                    if(inputReadIndex == 0): #read 1
                        currCluster.read1 = currRecord.seq
                        currCluster.qRead1 = currRecord.letter_annotations['phred_quality']
                    elif(inputReadIndex == 1): #read 2
                        currCluster.read2 = currRecord.seq
                        currCluster.qRead2 = currRecord.letter_annotations['phred_quality']
                    elif(inputReadIndex == 2): #index read 1
                        currCluster.index1 = currRecord.seq
                        currCluster.qIndex1 = currRecord.letter_annotations['phred_quality']
                    elif(inputReadIndex == 3): #index read2
                        currCluster.index2 = currRecord.seq
                        currCluster.qIndex2 = currRecord.letter_annotations['phred_quality']

                    #add/update the cluster
                    allClusterData[currRecord.id] = currCluster #add the new Cluster

                    #update the progress bar
                    update_progress(i,numRecordsInFile[inputReadIndex][inputFileIndex])
                    i = i + 1

                print 'len(allClusterData) = ' + str(len(allClusterData))

        #output data to file
        print 'writing output to file: ' + outputFilename
        with open(tempFilename(outputFilename),"w") as outputFileHandle: #temporary filename to indicate not finished
            lenAllClusterData = len(allClusterData)
            print 'writing ' + str(lenAllClusterData) + ' lines...'
            i = 0;
            for currID in allClusterData:
                cl = allClusterData[currID] #the current cluster
                #write to file
                outputFileHandle.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(currID, cl.filterID,
                                                                                              cl.read1, phredStr(cl.qRead1),
                                                                                              cl.read2, phredStr(cl.qRead2),
                                                                                              cl.index1,phredStr(cl.qIndex1),
                                                                                              cl.index2,phredStr(cl.qIndex2)))
                update_progress(i,lenAllClusterData)
                i = i + 1

        #rename ouput file to final filename indicating it is complete
            os.rename(tempFilename(outputFilename), outputFilename)

        #clean up database file
        if(useDatabase):
            os.remove(tempShelfFilename)

        if ('from_command_line' in keywordParams):
            #function was called from the command line
            return 0 #return exit value
        else:
            #function was imported and called from a library
            return outputFilename #return filename of output file

if __name__=='__main__':
    sys.exit(main(sys.argv, from_command_line=True))

