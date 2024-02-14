
import os
import sys
import time
import re
import uuid
import subprocess
import mergeFastqReadsToCPseq
import splitCPseqIntoTiles

from Bio import SeqIO

def mergeFastqs(fastqFolder, outputFilename):
    outArgs = ['mergeFastqReadsToCPseq']
    outArgs = outArgs + ['-rd',fastqFolder]
    outArgs = outArgs + ['-o',outputFilename]
    return mergeFastqReadsToCPseq.main(outArgs)

def splitCPseq(flowcellSide, CPseqFilename, outputFolder):
    outArgs = ['splitCPseqIntoTiles']
    outArgs = outArgs + ['-s',flowcellSide]
    outArgs = outArgs + ['-o',outputFolder]
    outArgs = outArgs + [CPseqFilename];
    return splitCPseqIntoTiles.main(outArgs)

def spawnMatlabJob(matlabFunctionCallString,globalVarsPath):
    try:
        #construct the command-line matlab call 
        functionCallString =                      "try,"
        functionCallString = functionCallString +     "addpath('{0}');".format(globalVarsPath) #placeholder TEMP DEBUG CHANGE
        functionCallString = functionCallString +     matlabFunctionCallString + ';'
        functionCallString = functionCallString + "catch e,"
        functionCallString = functionCallString +     "disp(getReport(e,'extended'));"
        functionCallString = functionCallString + "end,"
        functionCallString = functionCallString + "quit;"
    
        logFilename = 'matlabProcess_' + str(uuid.uuid4()) + str(time.time()) + '.tempLog' #timestamped logfile filename
    
        cmdString ='matlab -nodesktop -nosplash -singleCompThread -r "{0}"'.format(functionCallString)
        cmdString = cmdString + ' 1>> {0}'.format(logFilename)
        cmdString = cmdString + ' 2>> {0}'.format(logFilename)
       
        print 'issuing subprocess shell command: ' + cmdString
       
        returnCode = subprocess.call(cmdString,shell=True) #execute the command in the shell
        #returnCode2 = subprocess.call('stty sane',shell=True) #matlab messes up the terminal in a weird way--this fixes it 
    
        #read log file into a string
        try:
            with open(logFilename) as logFilehandle:
                logString = logFilehandle.read()
            # delete logfile
            try:
                os.unlink(logFilename)
            except OSError:
                pass
        except IOError:
            logString = 'Log file not generated for command "' + functionCallString + '".'
    
        # return log
        return logString
    except Exception, e:
        return 'Python exception generated in spawnMatlabJob: ' + e.message

def getMatlabCellArrayStringFromList(stringList):
    cellString = "{'"
    for currS in stringList[:-1]:
        cellString = cellString + currS + "','"
    cellString = cellString + stringList[-1] + "'}"
    return cellString

def alignmentFilter(CPseqFilename, filterFilename, outputFilename, globalVarsPath):
    try:
        matlabFunctionCallString = "alignmentFilterCPseq('{0}','{1}','{2}');".format(CPseqFilename,filterFilename,outputFilename)
        logString = spawnMatlabJob(matlabFunctionCallString,globalVarsPath)
        return (CPseqFilename,logString)
    except Exception,e:
        return(CPseqFilename,'Python excpetion generated in alignmentFilter: ' + e.message)
    
def getRegistrationOffset(CPseqFilename, imageFilename, dataScaling, filterSubsets, registrationOffsetMapFilename, maskImageFilename, globalVarsPath):
    try:
        filterSubsetString = getMatlabCellArrayStringFromList(filterSubsets)
        #matlabFunctionCallString = "GenerateRegistrationOffsetMap('{0}','{1}','{2}', {3}, '', '{4}','{5}');".format(CPseqFilename, imageFilename, dataScaling, filterSubsetString, registrationOffsetMapFilename, maskImageFilename)
        # V2 tries to use background subtraction for the registration offset map generation
        matlabFunctionCallString = "GenerateRegistrationOffsetMap_V2('{0}','{1}','{2}', {3}, '', '{4}','{5}');".format(CPseqFilename, imageFilename, dataScaling, filterSubsetString, registrationOffsetMapFilename, maskImageFilename)
        logString = spawnMatlabJob(matlabFunctionCallString,globalVarsPath)
        return (CPseqFilename,logString)
    except Exception,e:
        return(CPseqFilename,'Python excpetion generated in getRegistrationOffset: ' + e.message + e.stack)

def quantifyFluorescence(CPseqFilename, imageFilename, dataScaling, filterSubsets, registrationOffsetMapFilename, maskImageFilename, globalVarsPath,CPfluorFilename):
    try:
        filterSubsetString = getMatlabCellArrayStringFromList(filterSubsets)
        matlabFunctionCallString = "AnalyseImage('{0}','{1}','{2}', {3}, '{4}','{5}','{6}');".format(CPseqFilename, imageFilename, dataScaling, filterSubsetString, registrationOffsetMapFilename, maskImageFilename,CPfluorFilename)
        logString = spawnMatlabJob(matlabFunctionCallString,globalVarsPath)
        return (CPseqFilename,logString)
    except Exception,e:
        return(CPseqFilename,'Python excpetion generated in quantifyFluorescence: ' + e.message)

def getFilteredFilenameFromTileFilename(tileFilename,filteredDir):
    (path,filename) = os.path.split(tileFilename) #split the file into parts
    (basename,ext) = os.path.splitext(filename)
    return os.path.join(filteredDir,basename + '_filtered' + ext)

def getRoffFilenameFromImageFilename(tileImageFilename,roffDir):
    (path,filename) = os.path.split(tileImageFilename) #split the file into parts
    (basename,ext) = os.path.splitext(filename)
    return os.path.join(roffDir,basename + '.roff')

def getMaskFilenameFromImageFilename(tileImageFilename):
    (path,filename) = os.path.split(tileImageFilename) #split the file into parts
    (basename,ext) = os.path.splitext(filename)
    return os.path.join(path,basename + '.mask.tif')

def getQuantifiedFilenameFromCPseqFilename(CPseqFilename,CPfluor_dir):
    (path,filename) = os.path.split(CPseqFilename) #split the file into parts
    (basename,ext) = os.path.splitext(filename)
    return os.path.join(CPfluor_dir,basename + '.CPfluor')

def getQuantifiedFilenameFromCPseqFilenameAndImageStamp(CPseqFilename,tileImageFilename,CPfluor_dir):
    (imagePath,imageFilename) = os.path.split(tileImageFilename) #split tile image file into parts
    (imageBasename,imageExt) = os.path.splitext(imageFilename)
    stamp = imageBasename.split('_')[-1] #extract the stamp of the tile image (the text between the last underscore in the filename and the file extension)
    (path,filename) = os.path.split(CPseqFilename) #split the CPseq file into parts
    (basename,ext) = os.path.splitext(filename)
    return os.path.join(CPfluor_dir,basename + '_' + stamp + '.CPfluor')

def filenameMatchesAListOfExtensions(filename, extensionList):
    for currExt in extensionList:
        if filename.lower().endswith(currExt.lower()):
            return True
    return False

def getTileNumberFromFilename(inFilename):
    (path,filename) = os.path.split(inFilename) #split the file into parts
    (basename,ext) = os.path.splitext(filename)
    matches = re.findall('tile[0-9]{1,3}',basename.lower())
    tileNumber = ''
    if matches != []:
        tileNumber = '{:03}'.format(int(matches[-1][4:]))
    return tileNumber

def findTileFilesInDirectory(dirPath, extensionList, excludedExtensionList):
    dirList = os.listdir(dirPath)

    filenameDict = {}
    for currFilename in dirList:
        if filenameMatchesAListOfExtensions(currFilename,extensionList) and not filenameMatchesAListOfExtensions(currFilename,excludedExtensionList):
            currTileNum = getTileNumberFromFilename(currFilename)
            if currTileNum != '':
                filenameDict[currTileNum] = os.path.join(dirPath,currFilename)
    if len(filenameDict) == 0:
        print '      NONE FOUND'
    else:
        for currTile,currFilename in filenameDict.items():
            print '      found tile ' + currTile + ': "' + currFilename + '"'
    return filenameDict

def parseClusterID(clusterID):
    tokens = clusterID.split(':')
    
    flowcellID = ''
    tileID = ''
    if len(tokens) == 7:
        flowcellID = tokens[2] #get the full flowcell ID (3rd field of the cluster ID)
        tileID = tokens[4] #get the tile ID (5th field of the cluster ID)

    tokens = flowcellID.split('-')
    
    flowcellSN = ''
    if len(tokens) == 2:
        flowcellSN = tokens[1] #serial number is the second half after the dash of the flowcellID
        
    currSide = None
    currSwath = None  
    if (len(tileID) == 4) and stringIsValidInt(tileID): #for more modern clusterIDs with a side and swath (e.g. MiSeq data)
        currSide = int(tileID[0]) #which side of the flowcell (1= top surface, 2=bottom surface)
        currSwath = int(tileID[1]) #which swath on the surface

        #note - for the current  we will only use tiles from the top.  This code can be modified to output other configurations of tiles
        currTile = int(tileID[2:4]) #which tile number
    elif (len(tileID) >= 1) and (len(tileID) < 4) and stringIsValidInt(tileID): #for antiquated clusterIDs with just the tile number and no side and swath (e.g. GA data):
        currSide = 2
        currSwath = 1
        currTile = int(tileID)
    return (flowcellSN,currSide,currSwath,currTile)    

def getFirstFilenameInInputFilenames(inputFilenames):
    #go to the first vaild input filename
    firstFilename = ''
    for currFilenameList in inputFilenames:
        if currFilenameList != []:
            if currFilenameList[0] != []:
                firstFilename = currFilenameList[0]
    if not os.path.isfile(firstFilename):
        raise Exception('Input filename "' + firstFilename + '" is not a valid file!')
    return firstFilename

def getFlowcellSNfromClusterID(currID):
    (flowcellSN,currSide,currSwath,currTile) = parseClusterID(currID) #get the flowcell ID from the Cluster ID
    return flowcellSN

def getFlowcellSNfromFile(inFilename,filetype):
    try:
        currRecord = SeqIO.parse(inFilename, filetype).next() #read in one record
    except StopIteration:
        print 'No records were found in file "' + inFilename + '". Exiting...'
        sys.exit()
    flowcellSN = getFlowcellSNfromClusterID(currRecord.id)
    (path,filename) = os.path.split(inFilename) #split the (e.g. fastq) file into parts
    return (flowcellSN,path)

def findInputSequenceFilesInDirectory(read_dir,inputFilenames,**keywordParams):    
    #silent mode determines if files will be printed to the screen as they are found
    silentMode = False
    if ('silentMode' in keywordParams):
        silentMode = keywordParams['silentMode']
    
    if not os.path.isdir(read_dir):
        raise IOError('Input directory "' + read_dir + '" does not exist.')
    
    dirList = os.listdir(read_dir)
    numFound = 0
    r1Matches = [x for x in dirList if '_R1_' in x] # find files in the directory matching '_R1_'
    if r1Matches != []:
        if not silentMode:
            print 'read1 files found in directory "' + read_dir + '":'
        for currFilename in r1Matches:
            if currFilename not in inputFilenames[0]:
                inputFilenames[0].append(os.path.join(read_dir,currFilename))
                numFound = numFound + 1
                if not silentMode:
                    print '  ' + currFilename
    r2Matches = [x for x in dirList if '_R2_' in x] # find files in the directory matching '_R2_'
    if r2Matches != []:
        if not silentMode:
            print 'read2 files found in directory "' + read_dir + '":'
        for currFilename in r2Matches:
            if currFilename not in inputFilenames[1]:
                inputFilenames[1].append(os.path.join(read_dir,currFilename))
                numFound = numFound + 1
                if not silentMode:
                    print '  ' + currFilename
    i1Matches = [x for x in dirList if '_I1_' in x] # find files in the directory matching '_I1_'
    if i1Matches != []:
        if not silentMode:
            print 'index1 files found in directory "' + read_dir + '":'
        for currFilename in i1Matches:
            if currFilename not in inputFilenames[2]:
                inputFilenames[2].append(os.path.join(read_dir,currFilename))
                numFound = numFound + 1
                if not silentMode:
                    print '  ' + currFilename
    i2Matches = [x for x in dirList if '_I2_' in x] # find files in the directory matching '_I2_'
    if i2Matches != []:
        if not silentMode:
            print 'index2 files found in directory "' + read_dir + '":'
        for currFilename in i2Matches:
            if currFilename not in inputFilenames[3]:
                inputFilenames[3].append(os.path.join(read_dir,currFilename))
                numFound = numFound + 1
                if not silentMode:
                    print '  ' + currFilename
                    
    
    #if requiredToFind is True, one or more files must have been found in the directory
    requiredToFind = False
    if ('requiredToFind' in keywordParams):
        requiredToFind = keywordParams['requiredToFind']
    if requiredToFind:
        if numFound == 0: #if no input files were identified
            raise IOError("at least one input sequence (e.g. fastq) file must present in the read directory")
    return inputFilenames

def makeDirIfDoesntExist(dirName):
    if not os.path.exists(dirName):
        os.makedirs(dirName)


def splitPathIntoListOfDirectories(path):
    dirList = []
    while (path != '/') and (path != ''):
        (path,suffix) = os.path.split(path)
        dirList.append(suffix)
    return dirList[::-1] #return list in reverse order


#adds a '[<index>]' to the end of the filename if the file already exists to prevent overwriting
def renameIfFileExists(fullPath):
    if(os.path.exists(fullPath)):
        #file exists, rename
        (path,filename) = os.path.split(fullPath) #split the file into parts
        (basename,ext) = os.path.splitext(filename)
        matches = re.search('\[[0-9]*\]\Z',basename) #search to see if any indexes have already been added to the filename
        if matches == None:
            #if no indexes have already been added:
            outputFilename = os.path.join(path,basename + '[1]' + ext) #add a bracketed 1-index to the end of the filename
        else:
            index = int(basename[matches.start()+1:matches.end()-1])+1 # extract the index and increment it
            outputFilename = os.path.join(path,basename[0:matches.start()+1] + str(index) + ']' + ext) #assemble the new filename with the incremented index
        outputFilename = renameIfFileExists(outputFilename) #recurse to make sure the new filename we have created also doesn't already exist
        return outputFilename
    else:
        #file does not exist, make no changes
        return fullPath

def stringIsValidInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def tempFilename(fullPath): #appends a prefix to the filename to indicate the the file is incomplete
    (path,filename) = os.path.split(fullPath) #split the file into parts
    return os.path.join(path,'__' + filename)

def isTemporaryFilename(fullPath):
    (path,filename) = os.path.split(fullPath) #split the file into parts
    return filename.startswith('__')
    
def renamePathDict(filenameDict,newDir):
    for currKey,currFullPath in filenameDict.items():
        (currPath,currFilename) = os.path.split(currFullPath)
        filenameDict[currKey] = os.path.join(newDir,currFilename)
    return filenameDict
    
    
    
    
