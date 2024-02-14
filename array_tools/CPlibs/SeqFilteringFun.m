classdef SeqFilteringFun
    
    methods (Static)
        
        %------------------------------------------------------------------
        
        % Given a filter sequence, required-base bitstring, and a set of sequence
        % data, returns an array that has value 1 for every sequence that passed
        % the filter, and the value 0 for every sequence that failed the filter.

        % Peter McMahon
        % November 2013

        function passedFilter = filterSequences(readToFilter, filterSubsequence, filterBitString, sequences_read1, sequences_read2, sequences_indexread1, sequences_indexread2, startIndex, endIndex)
            if readToFilter == 1
                if(strcmpi(endIndex,'end'))
                    endIndex = length(sequences_read1{1});
                end
                sub_sequences_read1 = cellfun(@(x) x(startIndex:endIndex), sequences_read1, 'UniformOutput', false);
                passedFilter = cellfun(@(s) SeqFilteringFun.filterSequenceWithSubsequence(s, filterSubsequence, filterBitString), sub_sequences_read1);
            elseif readToFilter == 2
                if(strcmpi(endIndex,'end'))
                    endIndex = length(sequences_read2{1});
                end
                sub_sequences_read2 = cellfun(@(x) x(startIndex:endIndex), sequences_read2, 'UniformOutput', false);
                passedFilter = cellfun(@(s) SeqFilteringFun.filterSequenceWithSubsequence(s, filterSubsequence, filterBitString), sub_sequences_read2);
            elseif readToFilter == 3
                if(strcmpi(endIndex,'end'))
                    endIndex = length(sequences_indexread1{1});
                end
                sub_sequences_indexread1 = cellfun(@(x) x(startIndex:endIndex), sequences_indexread1, 'UniformOutput', false);
                passedFilter = cellfun(@(s) SeqFilteringFun.filterSequenceWithSubsequence(s, filterSubsequence, filterBitString), sub_sequences_indexread1);
            elseif readToFilter == 4
                if(strcmpi(endIndex,'end'))
                    endIndex = length(sequences_indexread2{1});
                end
                sub_sequences_indexread2 = cellfun(@(x) x(startIndex:endIndex), sequences_indexread2, 'UniformOutput', false);
                passedFilter = cellfun(@(s) SeqFilteringFun.filterSequenceWithSubsequence(s, filterSubsequence, filterBitString), sub_sequences_indexread2);
            else
                error('filterSequences:filterSequences', 'Invalid read index.');
            end
        end
        
        %------------------------------------------------------------------
        
        % Tests to see if a given sequence passes a filter: the filter is that a
        % given subsequence should match suitably well to the sequence, and
        % specific bases within the subsequence should exactly match the sequence.

        % Curtis Layton
        % November 2013
        % Updated: November 2013, PLM -- use 'glocal' option in NW, so that opening
        %                                and closing gap penalties are ignored;
        %                                updated to use "filter" nomenclature

        % Example:
        %   Input:
        %   sequence          = 'ATGCTGCCATGCATAAACGTACAT';
        %   filterSubsequence = 'ATGCATAAACGTAG';
        %   filterBitString   = '11110000000000'; % requires that 'ATGC' matches exactly; mismatches allowed for the rest of the sequence
        %
        %   [passedFilter, pValue, globalAlignment] = filterSequenceWithSubsequence(sequence, filterSubsequence, filterBitString)
        %
        %   Output:
        %   Since the following alignment is found:
        %
        %      ATGCTGCCATGCATAAACGTACAT
        %              |||||||||||||x  
        %      --------ATGCATAAACGTAG--
        %
        %   this sequence passes the filter. Note that as required, the first four
        %   bases of the filter subsequence match exactly.

        function [passedFilter, pValue, globalAlignment] = filterSequenceWithSubsequence(sequence, filterSubsequence, filterBitString)
            %check that filterSequence and filterBitString are the same length
            if(~isequal(size(filterSubsequence),size(filterBitString)))
                error('filterSequenceWithSubsequence:filterSequenceWithSubsequence','filter sequence and filter bitstring must be the same length');
            end
            if(~all((filterBitString == '1') | (filterBitString == '0')))
                error('filterSequenceWithSubsequence:filterSequenceWithSubsequence','filter bitstring must be all ones or zeros, representing required or unrequired bases, respectively')
            end
            %check that the sequences are not empty
            if (length(sequence) < 1) || (length(filterSubsequence) < 1)
                passedFilter = false;
                pValue = 1;
                globalAlignment = ['-';' ';'-']; %blank alignment
                return
            end
            
            % the sequence is assumed to pass the filter until it is found to fail
            passedFilter = true;

            %convert bitstring to logical array of bits
            filterBits = (filterBitString == '1'); 

            %#### STEP 1: initial local alignment to find matches by sequence similarity
            [pValue,globalAlignment] = SeqFilteringFun.p_nwalign(sequence,filterSubsequence); %global needleman-wunch alignment
            if(pValue > GlobalVars.pThreshold) %compare p-value to threshold 
                passedFilter = false; %if above threshold, it is not a match
                return;
            end     
            queryAlignString = globalAlignment(1,:); % get the alignment string for the query sequence
            filterAlignString = globalAlignment(3,:); % get the alignment string for the filter subsequence



            %#### STEP 2: check to see if gaps occur inside required regions (i.e. between two 1's in the bitstring

            filterGapsCollapsed = regexprep(filterAlignString,'-*','-'); %replace multiple dashes with a single, placeholder dash indicating the gap position, but discarding information about its length

            %find the indices of all positions before and after a gap, indexed in terms of filterBitstring and filterSubsequence
            filterRawIdx = find(filterGapsCollapsed == '-'); %find all position indices where gaps occur
            if(~isempty(filterRawIdx)) %if any gaps were found in the filter subsequence
                beforeIdxDiff = 1:length(filterRawIdx); %offset to relate gap positions to the index of filterBitstring before the gap
                afterIdxDiff = 0:length(filterRawIdx)-1; %offset to relate gap positions to the index of filterBitstring after the gap
                filterIdxBeforeGap = filterRawIdx - beforeIdxDiff; %calculate pre-gap indices for the filterBitstring
                filterIdxAfterGap = filterRawIdx - afterIdxDiff; %calculate post-gap indices for the filterBitstring
                if(filterIdxBeforeGap(1)==0) %if the first position is a gap we cannot check if the position prior is required, so discard it
                    if(length(filterIdxBeforeGap) > 1)
                        filterIdxBeforeGap = filterIdxBeforeGap(2:end);
                        filterIdxAfterGap = filterIdxAfterGap(2:end) - 1;
                    else
                        filterIdxBeforeGap = [];
                        filterIdxAfterGap = [];
                    end
                end
                if(~isempty(filterIdxAfterGap) && (filterIdxAfterGap(end)>length(filterBits))) %if the last position is a gap we cannot check if the position after is required, so discard it
                    if(length(filterIdxAfterGap) > 1)
                        filterIdxBeforeGap = filterIdxBeforeGap(1:end-1);
                        filterIdxAfterGap = filterIdxAfterGap(1:end-1);
                    else
                        filterIdxBeforeGap = [];
                        filterIdxAfterGap = [];
                    end
                end
                gapsInRequiredRegions = any((filterBits(filterIdxBeforeGap) & filterBits(filterIdxAfterGap))); %true if any position has a gap such that the position before and after are both required bases on the bitstring
                if(gapsInRequiredRegions)
                    passedFilter = false; %if gaps are between required bases on the filterSubsequence, it is not a match
                    return;
                end
            end



            %#### STEP 3: check to see if any substitutions occur in required regions

            %find the positions of the required bases, indexed in terms of the global alignment
            filterAlignGapPositions = (filterAlignString == '-'); %find all positions (in the uncollapsed alignment) where a gap occurs
            filterAlignIdx = zeros(1,length(filterSubsequence)); %initialize a vector to hold the indices to relate the filterBitstring positions to the alignment positions 
            cumulativeOffset = 0;
            origIdx = 1;
            for alignIdx = 1:length(filterAlignGapPositions)
                if(filterAlignGapPositions(alignIdx)) %if it is a gap position
                    cumulativeOffset = cumulativeOffset + 1; %increase the offset
                else
                    filterAlignIdx(origIdx) = origIdx + cumulativeOffset; %otherwise assign the offset index
                    origIdx = origIdx + 1;
                end
            end

            requiredIndices = find(filterBits); %find the required positions, indexed in terms of the filterSubsequence

            %check that all required positions are equal between the filter subsequence and the positions globally aligned to it, relating the filter subsequence positions to the alignment with filterAlignIdx
            subsInRequiredRegions = ~isequal(queryAlignString(filterAlignIdx(requiredIndices)),filterAlignString(filterAlignIdx(requiredIndices)));
            if(subsInRequiredRegions)
                passedFilter = false; %if any of the required bases aligned to the querySequence do not match the filterSubsequence, it is not a match
                return;
            end
        end
        
        %------------------------------------------------------------------
        
        %wrapper for a Needleman-Wunsch alignment with swalign, but returns a p-value that
        %is comparable between sequences of different lengths instead of a
        %highly-length-dependant and arbitrary alignment score

        %this p-value should be statistically valid for the comparison of two short
        %DNA sequences ~30nt using the Nuc44 scoring matrix (the default for
        %swalign when 'Alphabet' is 'NT').  Note that k and lambda must be re-fit
        %for any other other scoring matrices or scoring paramter sets besides
        %this.

        %for more info on the p-value caluclation from alignment score and how k and lambda were fit, see fitKandLambda.m and:
        %Stephen F. Altschul, Warren Gish, [27] Local alignment statistics, In: Russell F. Doolittle, Editor(s), Methods in Enzymology, Academic Press, 1996, Volume 266, Pages 460-480

        function  [pValue, alignment, startat, matrices, nwScore] = p_nwalign(querySequence, refSequence)
            % Needleman-Wunsch alignment 
            [nwScore, alignment, startat, matrices] = nwalign(querySequence,refSequence,'Alphabet','NT','Glocal', true);

            pValue = SeqFilteringFun.nwScoreToPvalue(nwScore, length(querySequence), length(refSequence), GlobalVars.EVD_k, GlobalVars.EVD_lambda);
        end

        % p value is the probabiity that this score or one higher occured by chance
        % for 2 sequences of lengths m and n
        function pValue = nwScoreToPvalue(nwScore, m, n, k, lambda)
            %m = effective length of first sequence
            %n = effective lenght of second sequence

            %mean
            u = log(k*m*n)/lambda;
            %p-value from extreme value distribution
            pValue = 1 - exp(-1*exp(-1*lambda*(nwScore-u)));
        end
        
        %------------------------------------------------------------------
        
        %generate many random DNA sequences and fit the resulting distribution of
        %alignment scores to obtain the parameters k and lambda that relate the
        %alignment score to the probability of obtaining that score from an
        %alignment of sequences of given lengths

        %for more info on these methods see:
        %Stephen F. Altschul, Warren Gish, [27] Local alignment statistics, In: Russell F. Doolittle, Editor(s), Methods in Enzymology, Academic Press, 1996, Volume 266, Pages 460-480

        function [k, lambda] = fitKandLambda()
            %parameters
            seqLength = 40; %length of the random sequences that will be generated
            numTrials = 1e5; %number of random sequences to generate
            nBins = 150; %number of bins for doing the score histogram

            scores = generateManyRandomAlignments(numTrials, seqLength); %generate many random alignments with associated scores
            [heights,centers] = hist(scores,nBins); %bin by score
            pdf = heights/numTrials; %normalize such that the histogram is a probability distribution

            %calculate the cumulative distribution
            cdf = pdf;
            for i = 2:nBins
                cdf(i) = cdf(i-1) + cdf(i);
            end

            %fit the cumulative distribution to an extreme value distribuion with
            %the parameters k and lambda
            fitParams = fitAlignmentScoreParams(seqLength, seqLength, centers, cdf);
            k = fitParams(1);
            lambda = fitParams(2);

            %plot the fit
            figure;
            plot(centers,cdf,'.','color','black');
            hold on;
            fitData = extremeValueCDF(fitParams,centers,seqLength,seqLength);
            plot(centers,fitData);
        end

        function fitParams = fitAlignmentScoreParams(m, n, x, cdf)

            %bounds and initial guesses
                           %     k    lambda
            lb           = [     0,      0];
            initialGuess = [  0.02,    0.2];
            ub           = [     1,     10];

            %optimization options
            options = optimset( 'TolX', 1e-18,... % TolX = convergence criteron on the change in parameters one iteration to the next
                            'TolFun', 1e-18,... % TolX = convergence criteron on the change in the value of the objective function one iteration to the next
                            'MaxFunEvals', 500,... %max number of times objective function can be evaluated
                            'Display', 'iter-detailed',... %display no output from optimization
                            'jacobian', 'off',... %use the jacobian
                            'DerivativeCheck', 'off',... %do not compare user-supplied jacobian (gradient of objective) to calculated finite-differencing derivatives
                            'LargeScale', 'on');

            %objective function is the extreme value cumulative distribution
            fun = @(params,x)extremeValueCDF(params,x,m,n);

            %fit
            fitParams = lsqcurvefit(fun,initialGuess,x,cdf,lb,ub,options);
        end

        %extreme value cumulative distribution
        function f = extremeValueCDF(params, x, m, n)
            %retreive params
            k = params(1);
            lambda = params(2);
            %mean
            u = log(k*m*n)/lambda;
            %function value
            f = exp(-1*exp(-1*lambda.*(x-u)));
        end

        %align many paris of randomly generated DNA sequences with and report scores
        function scores = generateManyRandomAlignments(numTrials, seqLength)
            scores = zeros(1,numTrials);

            for i = 1:numTrials
                sequence1 = generateRandomDNA(seqLength);
                sequence2 = generateRandomDNA(seqLength);
                scores(i) = nwalign(sequence1,sequence2,'Alphabet','NT');
            end
        end

        %generate random DNA sequences of a specified length
        function sequence = generateRandomDNA(length)
            bases = 'ATCG';
            sequence =  bases(randi(4,1,length));
        end
        
        %------------------------------------------------------------------
        
    end
end