% Load sequencing data from a MiSeq or GA run and output the distribution
% (histogram) of scores
%
% Input:  filename
%         format (can only be "CPseq")
%         readToFilter -- files contain up to four reads; this argument defines which one the procedure tries to align against (1 -- read1, 2 -- read2, 3 -- index read 1, 4 -- index read 2)
%         filterSequence (subsequence to check for in sequence)
%         num_bins (number of bins to use in the histogram)

% Output: histogram of NW scores for the given filter sequence against all
%         the sequences in the input data file
%
% Peter McMahon (pmcmahon@stanford.edu)
% October 2013
% November 2013 (updated for CPseq format)

function SeqDataNWDistribution(filename, format, readToFilter, filterSequence, num_bins)
    if strcmpi(format, 'CPseq')
        % file is in CPseq format
        num_lines = FileFun.countLinesInFile(filename);
        pValues = zeros(num_lines, 1); % allocate array that will store the NW p-value for each sequence (cluster)
        
        fin = fopen(filename);
        
        currentSeqIndex = 1;
        while (~feof(fin))
            [clusterID,filterID,read1,phredRead1,read2,phredRead2,indexRead1,phredIndexRead1,indexRead2,phredIndexRead2,num_sequences] = FileFun.readCPseqChunk(fin);
            
            if readToFilter == 1
                sequences = read1;
            elseif readToFilter == 2
                sequences = read2;
            elseif readToFilter == 3
                sequences = indexRead1;
            elseif readToFilter == 4
                sequences = indexRead2;
            else
                error('SeqDataNWDistribution:SeqDataNWDistribution', 'Invalid read to filter. Must be 1-4.');
            end
            
            currentChunk_pValues = cellfun(@(s) SeqFilteringFun.p_nwalign(s, filterSequence), sequences);
            
            pValues(currentSeqIndex:currentSeqIndex+num_sequences-1) = currentChunk_pValues;
            
            currentSeqIndex = currentSeqIndex + num_sequences;
        end
        
        fclose(fin);
    else
        error('SeqDataNWDistribution:SeqDataNWDistribution','Invalid file format. Must be "CPseq".');
    end
    
    % hist(log10(pValues), num_bins) % use MATLAB figure to display histogram
    [N,X] = hist(log10(pValues), num_bins);
    disp('Histogram of log10(pValues)');
    DisplayFun.displayHistogramASCII(N,X)
end