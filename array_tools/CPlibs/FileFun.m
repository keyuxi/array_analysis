classdef FileFun
    methods (Static)
        
        %loads an entire CPseq file into memory and puts into into a
        %ClusterData object
        function cDat = loadCPseq(filename, dataScaling)
            
            %count the number of lines in the file
%             fileLineCount = FileFun.countLinesInFile(filename);
            
            % open CPseq file for reading
            fid = fopen(filename);
            if(fid == -1)
                error('alignmentFilterCPseq:alignmentFilterCPseq',['could not open input file:' filename]);
            end
            disp(['reading from file ' filename ':']);
            
            %load data
%             disp(['loading ' num2str(fileLineCount) ' clusters...']);

            % parse the line on tabs
            %              1            2         3         4          5          6             7                 8                  9                  10
            % format: <cluster ID> <filterIDs> <read 1> <quality 1> <read 2> <quality 2> <index read 1> <quality index read 1> <index read 2> <quality index read 2>
            
            data = textscan(fid, '%s %s %s %s %s %s %s %s %s %s', 'Delimiter', '\t');

            clusterIDs = data{1};
            filterMatches = data{2};
%             read1 = data{3};
%             phredRead1 = data{4};
%             read2 = data{5};
%             predRead2 = data{6};
%             indexRead1 = data{7};
%             phredIndexRead1 = data{8};
%             indexRead2 = data{9};
%             phredIndexRead2 = data{10};

            % then parse the cluster id on colons
            %                            1            2            3          4        5         6        7
            % format of cluster ID: <machine id>:<run index>:<flowcell id>:<lane #>:<tile #>:<x coord>:<y coord>
            numLines = size(data{1},1);
            x = zeros(numLines,1);
            y = zeros(numLines,1);
            
            %parse the clusterID to get the x,y coords
            for i = 1:numLines
                data = textscan(clusterIDs{i},'%s %d %s %d %d %d %d', 'Delimiter', ':');
                x(i) = data{6};
                y(i) = data{7};
            end
            
            %scale data to match data source and images
            if(isKey(GlobalVars.scalingFactor,dataScaling))
                x = x / GlobalVars.scalingFactor(dataScaling);
                y = y / GlobalVars.scalingFactor(dataScaling);
            else
                errMsg = ['"' num2str(dataScaling) '" is not a valid data scaling factor. Valid choices are:   '];
                currKeys = keys(GlobalVars.scalingFactor);
                for i = 1:length(currKeys)
                    errMsg = [errMsg '"' currKeys{i} '"   '];
                end
                error('FileFun:loadCPseq',errMsg);
            end
            
            %rotate data if angle is given in GlobalVars file
            if(sum(ismember(properties(GlobalVars),'cluster_data_rotation_angle'))==1);
                angle = GlobalVars.cluster_data_rotation_angle;
                fprintf('Rotating data by %4.2f degrees\n', angle)
                
                x_unrotated = x;
                y_unrotated = y;
                
                % rotate x and y
                
                rot_matrix = [x_unrotated, y_unrotated]*[cosd(angle), -sind(angle); sind(angle), cosd(angle)];
                x = rot_matrix(:, 1);
                y = rot_matrix(:, 2);
            end
                

            filterMap = containers.Map(); %a map to hold all the filter names
            
            %parse the filter matches
            for i = 1:numLines
                filterMatches{i} = filterMatches{i}(filterMatches{i} ~= ' '); %remove any spaces
                if(~isempty(filterMatches{i}))
                    [cellArrayOfTokens, numTokens] = StringFun.tokenizeString(filterMatches{i}, ':');
                    for j = 1:numTokens
                        if ~isKey(filterMap,cellArrayOfTokens{j})
                            filterMap(cellArrayOfTokens{j}) = true;
                        end
                    end
                end
            end
            
            %create a ClusterData object with the data that was read in, and return
            cDat = ClusterData(clusterIDs, filterMatches, filterMap, x, y);
            
            %close CPfile
            fclose(fid);
        end

        function [clusterID,filterID,read1,phredRead1,read2,phredRead2,indexRead1,phredIndexRead1,indexRead2,phredIndexRead2,num_sequences] = readCPseqChunk(fin_seqdata)
            if(fin_seqdata == -1)
                error('FileFun:readCPseqChunk','input file not open');
            end
            
            data = textscan(fin_seqdata, '%s %s %s %s %s %s %s %s %s %s', GlobalVars.chunkSize, 'Delimiter', '\t'); % read in chunkSize lines (at most), and parse the lines
            %              1            2         3         4          5          6             7                 8                  9                  10
            % format: <cluster ID> <filterIDs> <read 1> <quality 1> <read 2> <quality 2> <index read 1> <quality index read 1> <index read 2> <quality index read 2>
            num_sequences = size(data{1},1); % number of sequences in this chunk that we have just read (will be at most chunkSize, but could be fewer)

            clusterID = data{1};
            filterID = data{2};
            read1 = data{3};
            phredRead1 = data{4};
            read2 = data{5};
            phredRead2 = data{6};
            indexRead1 = data{7};
            phredIndexRead1 = data{8};
            indexRead2 = data{9};
            phredIndexRead2 = data{10};
        end
        
        function writeCPseqChunkToFile(fid, clusterID, filterID, read1, phredRead1, read2, phredRead2, indexRead1, phredIndexRead1, indexRead2, phredIndexRead2)
            %check length of data
            if(~isequal(size(clusterID,1),size(filterID,1),size(read1,1),size(phredRead1,1),size(read2,1),size(phredRead2,1),size(indexRead1,1),size(phredIndexRead1,1),size(indexRead2,1),size(phredIndexRead2,1)))
                error('alignmentFilterCPseq:writeChunkToFile','all data must be of the same length');
            end

            numLines = size(clusterID,1);

            %allocate a cell array to hold the output chunk
            outputChunk = cell(numLines,1);

            for i = 1:numLines
                outputChunk{i} = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', clusterID{i}, filterID{i}, read1{i}, phredRead1{i}, read2{i}, phredRead2{i}, indexRead1{i}, phredIndexRead1{i}, indexRead2{i}, phredIndexRead2{i});
            end

            %write chunk to file
            if(fid ~= -1)
                fprintf(fid, '%s', outputChunk{:}); %write the transpose to the file
            else
                error('FileFun:writeChunkToFile','file not open');
            end
        end

        function saveCPfluorFile(cDat, CPfluorFilename)
            %check file extension
            [CPfluorPath, CPfluorName, CPfluorExt] = fileparts(CPfluorFilename);
            if(~strcmp(lower(CPfluorExt),'.cpfluor'))
                CPfluorFilename = [CPfluorFilename '.CPfluor'];
            end

            %write data to the cell array that will be written to file
            outputChunk = cell(cDat.numClusters,1);
            for i = 1:cDat.numClusters
                %format: <Cluster ID>:<fit flag>:<amplitude>:<sigma>:<y>:<x>
                outputChunk{i} = sprintf('%s:%d:%f:%f:%f:%f\n', cDat.IDs{i},cDat.alreadyFit(i),cDat.fitAmplitudes(i),cDat.fitSigmas(i),cDat.fitY(i),cDat.fitX(i));
            end
            
            %write to file
            fid = fopen(CPfluorFilename,'wt'); %open file
            if(fid == -1)
                error('FileFun:writeChunkToFile','could not open file "%s" for output', CPfluorFilename);
            else
                fprintf(fid, '%s', outputChunk{:}); %write the data to the file
                fclose(fid); %close file
            end
        end
       
        function lineCount = countLinesInFile(filename)
            fh = fopen(filename, 'rt');
            assert(fh ~= -1, 'Could not read: %s', filename);
            x = onCleanup(@() fclose(fh));
            lineCount = 0;
            while ~feof(fh)
                lineCount = lineCount + sum( fread( fh, 16384, 'char' ) == char(10) );
            end
        end
        
	% Load colormap from .CPfiltercolormap file. Associates filter names with colors; when cluster data are plotted, 
	% points from sets matching different filters will be assigned colors based on the filter name. A special "ELSE"
	% label can be used to specify the color used to plot points from all filters that don't have specified colors,
	% or clusters that didn't match any filter.
	%
	% e.g. *** test.CPfiltercolormap ***
	%      FilterName1 red
	%      FilterName2 green
	%      ELSE blue
        function filterColormap = loadFilterColormap(filename)
            if isempty(filename)
                filterColormap = containers.Map();
            else
                fid = fopen(filename, 'rt');

                data = textscan(fid, '%s %s', 'Delimiter', '\t');
                filterNames = data{1};
                filterColors = data{2};

                filterColormap = containers.Map(filterNames, filterColors);

                if ~isKey(filterColormap, 'ELSE')
                    filterColormap('ELSE') = 'blue'; % if the CPfiltercolormap file did not specify a color for "ELSE", then define it as blue
                end

                fclose(fid);
            end
	end


    end % END static methods
end
