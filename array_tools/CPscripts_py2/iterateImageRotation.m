% originally written by Sarah Denny
%
% Goal: try different rotations of data and evaluate correlation peak success.
% Method: rotate image, register, and evaluate success of registration by correlation peak height.
% 
% Note: initial default setting is to try broadly space angles (+/- 3 degrees,
% 10 iterations. May want to input your own set after first iteration.

% modified by CJL 01.30.2015


function [angles, offsets_x, offsets_y,peak_heights] = iterateImageRotation(seqDataFilename, imageFilename, dataScaling, filterSubsets, angles)

    %% load image
    disp(['Loading image from file ' imageFilename ':'])
    rawImage = double(imread(imageFilename));
   
    %load sequence data
    disp('Loading and scaling sequence data...')
    cDat = FileFun.loadCPseq(seqDataFilename, dataScaling);
    
    %register image
    disp('Performing registration of the data to the image...')
    
    %mark the subset of clusters that will be used for registration
    cDat.currSubset = RegistrationFun.identifyFilterSubsetsForRegistration(cDat, filterSubsets);
    
    %define the central tile to register (size set in GlobalVars)
    centralSubtileCoords = ImageFun.getCentralSubtileCoords(size(rawImage));

    % change angles of rotation by small iterations if not given as input
    if(~exist('angles', 'var'))
        angles = -2:0.5:2; % in degrees
    end
    
    % initiate arrays to save offsets
    num_iter = length(angles);
    offsets_x = zeros(num_iter, 1);
    offsets_y = zeros(num_iter, 1);
    peak_heights = zeros(num_iter, 1);

    %define the bounds of the data to register
    dataBounds = [min(cDat.currX()) min(cDat.currY()) max(cDat.currX()) max(cDat.currY())];  
    
    %only use points in the subtile and that are part of the
    %current subset (as defined by filter matching) for registration
    subtileAdjY = cDat.adjY(cDat.currSubset);
    subtileAdjX = cDat.adjX(cDat.currSubset);
    
    for i=1:num_iter;

        % rotate data
        a = angles(i);
        rotImage = imrotate(rawImage, -a); %invert rotation so that the recorded angle will be the correct input angle into GlobalVars.cluster_data_rotation_angle

        % crop out a subtile region of the real image for correlation
        subtileImage = imcrop(rotImage, centralSubtileCoords);
        subtileOffset = [(1 - centralSubtileCoords(2)) (1 - centralSubtileCoords(1))]; %the offset of the subimage relative to the raw image
        
        % define the bounds of the synthetic image such that they are guaranteed
        % to be larger than the central region of the real image (normxcorr2 demands that
        % the first argument (real image) have smaller dimensions than the second argument(synth image))
        synthImageSizeY = max( ceil(max(subtileAdjY)-min(subtileAdjY))+1 , centralSubtileCoords(4)+2 );
        synthImageSizeX = max( ceil(max(subtileAdjX)-min(subtileAdjX))+1 , centralSubtileCoords(3)+2 );

        %create a synthetic image for cross-correlation
        synthImg = ImageFun.makeSynthImagePoints(subtileAdjY, subtileAdjX, synthImageSizeY, synthImageSizeX, GlobalVars.synthIntensityValue);
        
        % Normalized crosscorrelation of the subtile image versus the synthetic image
        corrMatrix = normxcorr2(subtileImage, synthImg);

        %fit the correlation matrix to a gaussian to find the peak
        %correlation to sub-pixel resolution.  The correlation peak height is
        %used as a measure of registration success
        [currTileRegistrationSuccess yPeak xPeak correlationPeakHeight] = RegistrationFun.fitCorrelationGaussian(corrMatrix, [Inf Inf]);
        
        % sum the relative offsets of the position of the central subimage (synthetic data) registered onto the image
        corrOffset = [(yPeak - size(subtileImage,1)) (xPeak - size(subtileImage,2))]; %the offset of the correlation peak relative to the subtile image size
        currOffsetY = min(subtileAdjY) + corrOffset(1) + subtileOffset(1) - 1;
        currOffsetX = min(subtileAdjX) + corrOffset(2) + subtileOffset(2) - 1;

        % register and record 
        offsets_x(i) = currOffsetX;
        offsets_y(i) = currOffsetY;
        peak_heights(i) = correlationPeakHeight;
    end

    %% plot offsets and peak heights
    figure()
    subplot(2,1,1)
    labels = cellstr( num2str(angles', '%4.4f') );
    plot(offsets_x, offsets_y, 'ro')
    text(offsets_x, offsets_y, labels, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')
    xlabel('offset in x')
    ylabel('offset in y')
    
    subplot(2,1,2)
    plot(angles, peak_heights, 'bo-')
    xlabel('rotation angle')
    ylabel('peak height')
end
