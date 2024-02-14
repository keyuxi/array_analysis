%check if a tile image registeres to tile data (.CPseq format), return the offset and a image of the registration

function  [success, offsetY, offsetX, centeredTilePositionOffsetY, centeredTilePositionOffsetX, outputImageFilename, zoomOutputImageFilename] = checkTileRegistration(seqDataFilename, imageFilename, dataScaling, filterSubsets, imageOutputPath, remoteCall)

    checkRegUUID = '9d7bd2bb-5f26-4893-89ac-641e30ec8ab9'; %unique string to identify correct output when called remotely

    %### DEBUG ###
    if(GlobalVars.debug)
        global debugFigMan;
        if(isprop(debugFigMan,'globallyDismissed'))
            debugFigMan.delete();
        end
        debugFigMan = debugFigureManager();
    end
    %### END DEBUG ###

    %the filter subsets that will be used for registration
    if(~exist('filterSubsets','var'))
        filterSubsets = []; %an empty filter subset name list means all clusters in the file will be used for registration
    end
    
    if(~exist('remoteCall','var'))
        remoteCall = false; %indicates if the program was called from the command line by a remote program -- only affects printed output formatting
    end
    
    %output QC image files
    generateImage = true;
    outputImageFilename = 'NONE'; %flag indicating the output image file was not generated
    zoomOutputImageFilename = 'NONE'; %flag indicating the output image file was not generated
    
    if ~exist('imageOutputPath','var') || isempty(imageOutputPath)
        imageOutputPath = []; %an empty output image path means that the registration QC images will not be generated
        generateImage = false;
    else
        if(~exist(imageOutputPath,'dir'))
            %output path for image output is invalid
            outputImageFilename = 'ERR'; %flag indicating there was an error generating the output file
            zoomOutputImageFilename = 'ERR'; %flag indicating there was an error generating the output file
            generateImage = false;
        end
    end

    %create colorMap for filterSubsets (to visualize the clusters being used for registration)
    if ~isempty(filterSubsets)
        filterColors = cell(length(filterSubsets),1);
        for idx = 1:length(filterSubsets)
            filterColors{idx} = 'green'; % assign each filter in the registration subset to be green
        end
        filterColormap = containers.Map(filterSubsets, filterColors);
        filterColormap('ELSE') = 'blue'; % assign filters not used in registration to be blue
    else
        filterColormap = containers.Map(); % if all clusters are to be used for registration make empty colorMap.
                                           % The plot function will assign colors automatically.
    end

    %read in image
    disp(['Loading image from file ' imageFilename ':'])
    rawImage = double(imread(imageFilename));
    
    % remove background from the image with morphological opening
    % (this greatly improves cross correlation performance)
    bgImage = ImageFun.getBackgroundFromClusterImage(rawImage);
    bgSubtractedImage = rawImage - bgImage;
    
    %get the center coordinates of the image
    imgSize = size(bgSubtractedImage);
    imageCenter = [(imgSize(1)/2) (imgSize(2)/2)];
    
    %load sequence data
    disp('loading and scaling sequence data...');
    cDat = FileFun.loadCPseq(seqDataFilename, dataScaling);
    cDat.currSubset = RegistrationFun.identifyFilterSubsetsForRegistration(cDat, filterSubsets); %mark the subset of clusters that will be used for registration
        
    %generate subtile boundaries -- same number of tiles as the coarse step of hierarchal registration
    subtileMap = ImageFun.getSubtileMapCoords(size(bgSubtractedImage), GlobalVars.nTilesRegistrationCheck(1), GlobalVars.nTilesRegistrationCheck(2));
       
    %register subtiles against the whole image (unbounded)
    success = false;
    offsetY = 0;
    offsetX = 0;
    centeredTilePositionOffsetY = 0;
    centeredTilePositionOffsetX = 0;
    offsetListY = [];
    offsetListX = [];

    %iterate over all subtiles
    for i = 1:size(subtileMap,1)
        %get the coords of the current subtile image
        subImageRect = subtileMap(i,:);
        
        % crop out a subtile region of the real image for correlation
        subtileImage = imcrop(rawImage, subImageRect + [GlobalVars.cropBorder GlobalVars.cropBorder -(2*GlobalVars.cropBorder) -(2*GlobalVars.cropBorder)]); %crop the border of the subtile further
        subtileOffset = [(1 - subImageRect(2) - GlobalVars.cropBorder) (1 - subImageRect(1) - GlobalVars.cropBorder)]; %the offset of the subimage relative to the raw image

        %only use points in the subtile and that are part of the
        %current subset (as defined by filter matching) for registration
        subtileAdjY = cDat.adjY(cDat.currSubset);
        subtileAdjX = cDat.adjX(cDat.currSubset);

        % define the bounds of the synthetic image such that they are guaranteed
        % to be larger than the central region of the real image (normxcorr2 demands that
        % the first argument (real image) have smaller dimensions than the second argument(synth image))
        synthImageSizeY = max( ceil(max(subtileAdjY)-min(subtileAdjY))+1 , subImageRect(4)+2 );
        synthImageSizeX = max( ceil(max(subtileAdjX)-min(subtileAdjX))+1 , subImageRect(3)+2 );

        %create a synthetic image for cross-correlation
        synthImg = ImageFun.makeSynthImagePoints(subtileAdjY, subtileAdjX, synthImageSizeY, synthImageSizeX, GlobalVars.synthIntensityValue);
        
        % Normalized crosscorrelation of the subtile image versus the synthetic image
        corrMatrix = normxcorr2(subtileImage, synthImg);

        %fit the correlation matrix to a gaussian to find the peak
        %correlation to sub-pixel resolution.  The correlation peak height is
        %used as a measure of registration success
        [currTileRegistrationSuccess yPeak xPeak correlationPeakHeight] = RegistrationFun.fitCorrelationGaussian(corrMatrix, [Inf Inf]);
        
        % sum the relative offsets of the position of the central subimage (synthetic data) registered onto the actual image
        corrOffset = [(yPeak - size(subtileImage,1)) (xPeak - size(subtileImage,2))]; %the offset of the correlation peak relative to the subtile image size
        currOffsetY = min(subtileAdjY) + corrOffset(1) + subtileOffset(1) - 1;
        currOffsetX = min(subtileAdjX) + corrOffset(2) + subtileOffset(2) - 1;
                    
        %success if multiple subtiles register at a threshold in the same frame
        if correlationPeakHeight > GlobalVars.corrSuccessPeakHeightRegistrationCheck
        
            %add registration offsets that pass threshold to a list for cross-comparison
            offsetListY(length(offsetListY)+1) = currOffsetY;
            offsetListX(length(offsetListX)+1) = currOffsetX;

            %if more than two subtiles register above threshold to the same place, we have found the correct registration
            withinToleranceList = getWithinTolerance(currOffsetY, currOffsetX, offsetListY, offsetListX, GlobalVars.matchingCorrelationTolerance);

            if length(withinToleranceList) >= GlobalVars.numTileToMatch
                offsetY = mean(offsetListY(withinToleranceList));
                offsetX = mean(offsetListX(withinToleranceList));
                success = true;
                break;
            end
        end
    end
    
    if success
        %move the cluster positions
        cDat.adjY = cDat.adjY - offsetY;
        cDat.adjX = cDat.adjX - offsetX;

        %get the center coordinates of the data
        minDataY = min(cDat.adjY(:));
        maxDataY = max(cDat.adjY(:));
        minDataX = min(cDat.adjX(:));
        maxDataX = max(cDat.adjX(:));
        dataCenter = [minDataY+((maxDataY - minDataY)/2) minDataX+((maxDataX - minDataX)/2)];
        
        %calcuate the centered tile position offset
        centeredTilePositionOffsetY = imageCenter(1)-dataCenter(1);
        centeredTilePositionOffsetX = imageCenter(2)-dataCenter(2);
        
    end
    
    %### DEBUG ###
    if(GlobalVars.debug)
        showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:RegisterHierarchical: Would you like to see the global registration?');
        if(showDebug)
              PlotFun.plotGlobalRegistrationByFilter(rawImage, cDat, [], [], 'global registration', filterColormap);
        end
        drawnow;
        debugFigMan.clearLocal(); %clear debug dismissals from the previous level of registration
    end
    %### END DEBUG ###
    
    %plot registration for visual confirmation and save as images
    if success && generateImage

        %create filenames for the output images
        [pathstr,name,ext] = fileparts(imageFilename);
        outputImageFilename = fullfile(imageOutputPath, [name '_Reg.tif']);
        zoomOutputImageFilename = fullfile(imageOutputPath, [name '_RegZoom.tif']);
        
        figureHandle = figure('visible','off'); %create a figure, but do not display it in the gui
        axesHandle = axes('Parent',figureHandle);

        %global position,  zoomed out
        currTitle = sprintf('Image Tile Registration to %s',seqDataFilename);
        PlotFun.plotGlobalRegistrationByFilter(rawImage, cDat, [], [], currTitle, filterColormap, figureHandle, axesHandle); %plot the registration in the hidden figure
        print(figureHandle, '-r80', '-dtiff', outputImageFilename);  %save the figure as an image
        
        %zoomed in
        currTitle = sprintf('Zoomed Image Tile Registration to %s',seqDataFilename);
        title(axesHandle, currTitle);
        zoom(figureHandle,8); %zoom in
        print(figureHandle, '-r80', '-dtiff', zoomOutputImageFilename);  %save the figure as an image
        
        close(figureHandle); %close figure
    end
    
    %if called remotely, after a unique identifier (for parsing the output), output space-delimited results 
    if remoteCall
        finalOutput = sprintf('%s %d %f %f %f %f %s %s\n', checkRegUUID, success, offsetY, offsetX, centeredTilePositionOffsetY, centeredTilePositionOffsetX, outputImageFilename, zoomOutputImageFilename);
        disp(finalOutput);
    else
        %output formatted for human readibility
        disp('')
        if success
            disp('Registration FOUND');
            msg = sprintf('   Y offset = %f', offsetY); disp(msg);
            msg = sprintf('   X offset = %f', offsetX); disp(msg);
            msg = sprintf('   offset in pixels from current position to centered tile position (Y,X) = (%f, %f)', centeredTilePositionOffsetY, centeredTilePositionOffsetX); disp(msg);
            msg = sprintf('   global view registration image = %s', outputImageFilename); disp(msg);
            msg = sprintf('   zoomed-in registration image = %s', zoomOutputImageFilename); disp(msg);
        else
            disp('Registration not found')
        end
    end
end

%test if a point is within tolerance of any other point in the list
function withinTolerance = getWithinTolerance(currOffsetY, currOffsetX, offsetListY, offsetListX, tolerance)
    n = min(length(offsetListY),length(offsetListX));
    withinTolerance = [];
    for i = 1:n
        if (abs(currOffsetY - offsetListY(i)) < tolerance) && (abs(currOffsetX - offsetListX(i)) < tolerance)
            withinTolerance(length(withinTolerance)+1) = i;
        end
    end
end