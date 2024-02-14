% Load a .CPseq sequencing data file, filter it by sequence alignment to reference sequences using a set of
% filters, and output all the data into a new .CPseq file, adding the filter IDs.
%
% Input:  input_filename
%         inputFilenameFilterExpressions (filename of file containing filtering expressions)
%
% Output: files in same format is the input (.CPseq), but with the names of matching filters added to the <filter IDs> field.
%
% Example usage:
%   > alignmentFilterCPseq('test.CPseq','filter_seq_example.CPfilter')
%
% Curtis Layton (curtis.layton@stanford.edu) and Peter McMahon (pmcmahon@stanford.edu)
% November 2013


% .CPseq format (tab delimited):
%
%              1            2         3         4          5          6             7                 8                  9                  10
% format: <cluster ID> <filterIDs> <read 1> <quality 1> <read 2> <quality 2> <index read 1> <quality index read 1> <index read 2> <quality index read 2>
%
%                            1            2            3          4        5         6        7
% format of cluster ID: <machine id>:<run index>:<flowcell id>:<lane #>:<tile #>:<x coord>:<y coord>


function alignmentFilterCPseq(inputFilenameSeqdata, inputFilenameFilterExpressions, outputFilename)
%     profile on;
    
    %assemble the output filename, if one was not passed in
    if(~exist('outputFilename','var'))
        [seqPath, seqName, seqExt] = fileparts(inputFilenameSeqdata); %dissect the input filename
        outputPath = seqPath; %default output path is the same as the input
        if ~isempty(outputPath) && (outputPath(end) ~= StringFun.pathSlash()) %add a trailing slash to the output path if it does not already exist
            outputPath = [outputPath StringFun.pathSlash()];
        end
        outputFilename = [outputPath seqName '_filtered.CPseq']; %concatenate a "_filtered" designation onto the end of the input filename
    end
    
    %count the number of lines in the file
    inputFileLineCount = FileFun.countLinesInFile(inputFilenameSeqdata);
    
    % load filter sequences and bitstrings
    fin_filterSequences = fopen(inputFilenameFilterExpressions);
    filterFileLines = textscan(fin_filterSequences, '%s', 'Delimiter', '\n'); % one filter expression per line
    filterFileLines = filterFileLines{1};
    numFilters = length(filterFileLines); % number of sequences that we should filter the sequence data for
    filterNames = cell(1,numFilters);
    filterExpressions = cell(1,numFilters);
    
    %extract filter names and logical expressions from each line of the filter file
    for i = 1:numFilters
        [lineTokens, numTokens] = StringFun.tokenizeString(filterFileLines{i}, ':');
        if(numTokens ~= 2)
            error('alignmentFilterCPseq:alignmentFilterCPseq','Invalid filter file. Example filter file format: "<filter name>:filter(<read idx>,''<sequence>'',''<bitstring>'')". Filter names must not contain any colons.');
        end
        filterNames{i} = lineTokens{1};
        if(~isValidFilterName(filterNames{i}))
            errorMsg = ['invalid filter name "' filterNames{i} '".  Filter names must be alphanumeric and may contain underscores, dashes, and periods, but no other characters'];
            error('alignmentFilterCPseq:alignmentFilterCPseq',errorMsg);
        end
        filterExpressions{i} = lineTokens{2};
    end

    % sequencing data input file is in cp format
    fin_seqdata = fopen(inputFilenameSeqdata);
    if(fin_seqdata == -1)
        error('alignmentFilterCPseq:alignmentFilterCPseq',['could not open input file:' inputFilenameSeqdata]);
    end
    disp(['reading from file ' inputFilenameSeqdata ':']);
    
    % open the output file
    fout_seqdata = fopen(outputFilename,'wt');
    if(fout_seqdata == -1)
        error('alignmentFilterCPseq:alignmentFilterCPseq',['could not open output file:' outputFilename]);
    end
    disp(['writing to file ' outputFilename ':']);

    %init progress bar
    disp(['filtering ' num2str(inputFileLineCount) ' entries by sequence alignment...']);
    progbar = DisplayFun.progressBar(inputFileLineCount);
    
    %read through the sequence input file in chunks
    lineCount = 0;
    while (~feof(fin_seqdata))

        [currClusterID,currFilterID,currRead1,currPhredRead1,currRead2,currPhredRead2,currIndexRead1,currPhredIndexRead1,currIndexRead2,currPhredIndexRead2,num_sequences] = FileFun.readCPseqChunk(fin_seqdata);

        lineCount = lineCount + num_sequences;
          
        filter = @(readToFilter, filterSequence, filterBitString, startIndex, endIndex)SeqFilteringFun.filterSequences(readToFilter, filterSequence, filterBitString, currRead1, currRead2, currIndexRead1, currIndexRead2, startIndex, endIndex); % define an anonymous function that allows the filterSequences function to be called without needing to explicitly pass the sequences_* variables
        
        for filterFileIndex = 1:numFilters
            % for each filter, check if the sequences in the current chunk pass that filter

            currentFilterExpression = filterExpressions{filterFileIndex};

            eval(['passedFilter = ' currentFilterExpression ';']);

            filteredIndices = find(passedFilter);
            
            for i = 1:length(filteredIndices)
                if(isempty(currFilterID{filteredIndices(i)}))%if there is no filter ID
                    currFilterID{filteredIndices(i)} = filterNames{filterFileIndex}; %add the current filter ID
                else%if there is already an existing filter ID
                    tokens = StringFun.tokenizeString(filterNames{filterFileIndex},':'); %split the filterID string into individual filterIDs
                    if(isempty(find(strcmp(currFilterID{filteredIndices(i)},tokens), 1)))%make sure this seqeunce doesn't already have the passing filter name in its filter ID
                        currFilterID{filteredIndices(i)} = [currFilterID{filteredIndices(i)} ':' filterNames{filterFileIndex}]; %if not present, append the current filter ID
                    end
                end
            end
        end
        
        % Write data to output file, adding matching filters to the <filter IDs> field
        FileFun.writeCPseqChunkToFile(fout_seqdata, currClusterID, currFilterID, currRead1, currPhredRead1, currRead2, currPhredRead2, currIndexRead1, currPhredIndexRead1, currIndexRead2, currPhredIndexRead2);

        %update progress bar
        progbar(lineCount,round(100/GlobalVars.chunkSize)+1);
    end
    fprintf('\n'); %print newline for progress bar

    fclose(fin_seqdata);
    fclose(fout_seqdata);
     
%     profile viewer;
%     profile off;
end


%only allow a restricted character set for filter names.  These names will go into filenames.  Note the conspicuous omission of spaces and other whitespace
function valid = isValidFilterName(name)
    valid = all((name>='a'&name<='z') | (name>='A'&name<='Z') | (name>='0'&name<='9') | (name=='_') | (name=='-') | (name=='.'));
end