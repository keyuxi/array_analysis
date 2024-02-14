%%% Uses lsqcurvefit to fit Gaussians over sub-images  %%%
function [cDat,bgImage] = QuantifyImage(I,cDat,zeroPix)
    %### DEBUG ###
    if(GlobalVars.debug)
        global debugFigMan;
    end
    profile on;
    %### END DEBUG ###
    
    %we will tile over the image I and fit 2D gaussians to each data point 
    %that falls within each small subtile.  The subtiles overlap each other by 
    %(bound) pixels to prevent boundary effects

    %fit all unexcised data
    unexcisedIndices = find(~cDat.excised);
    adjY = cDat.adjY(unexcisedIndices);
    adjX = cDat.adjX(unexcisedIndices);
    
    %define function and set optimization parameters
    fun = @FitFun.multiGaussian_2D; % objective function for optimization (see below)

    options = optimset( 'TolX', GlobalVars.multiGauss_TolX,... % TolX = convergence criteron on the change in parameters one iteration to the next
                        'TolFun', GlobalVars.multiGauss_TolFun,... % TolFun = convergence criteron on the change in the value of the objective function one iteration to the next
                        'MaxFunEvals', GlobalVars.multiGauss_MaxFunEvals,... %max number of times objective function can be evaluated
                        'Display', 'off',... %display no output from optimization
                        'jacobian', 'on',... %use the jacobian
                        'DerivativeCheck', 'off',... %do not compare user-supplied jacobian (gradient of objective) to calculated finite-differencing derivatives
                        'PrecondBandWidth', 0,...
                        'LargeScale', 'on');

    %round position vectors outside the loop to avoid repetitious rounding inside the loop
    roundedAdjY = round(adjY);
    roundedAdjX = round(adjX);

    bgImage = zeros(size(I)); %initialize the background image (which will be filled in with the background fit levels)

    %define the coords of all the subtiles that will be fit
    [innerSubtileMap outerSubtileMap, numSubtiles] = ImageFun.getOverlappingSubtileMapCoords(size(I), GlobalVars.innerWindowSize, GlobalVars.innerWindowSize, GlobalVars.overlap, GlobalVars.overlap);
    
    % Make x,y coordinates for every pixel in the sub region. meshgrid creates 2D coordinate grid arrays based on the size and value of the input vectors
    maxOuterSizeX = max(outerSubtileMap(:,3));
    maxOuterSizeY = max(outerSubtileMap(:,4));
    [fitWindowX,fitWindowY] = meshgrid(1:maxOuterSizeX,1:maxOuterSizeY);  

    %init progress bar
    disp('Fitting...');
    progbar = DisplayFun.progressBar(numSubtiles);

    %loop through overlapping subtiles and fit
    for subtileIdx = 1:numSubtiles
                
        %define the dimensions of the current overlapped subtile 
        %inner subtile (fits from clusters within this region will be kept)
        currInnerSubtileMap = innerSubtileMap(subtileIdx,:);
        innerMinX = currInnerSubtileMap(1);
        innerMinY = currInnerSubtileMap(2);
        innerWidthX = currInnerSubtileMap(3);
        innerHeightY = currInnerSubtileMap(4);
        innerMaxX = innerMinX + innerWidthX - 1;
        innerMaxY = innerMinY + innerHeightY - 1;
        %outer subtile (the boundary region is included in the fit to preven edge effects, but the fits from clusters in this region are discarded)
        currOuterSubtileMap = outerSubtileMap(subtileIdx,:);
        outerMinX = currOuterSubtileMap(1);
        outerMinY = currOuterSubtileMap(2);
        outerWidthX = currOuterSubtileMap(3);
        outerHeightY = currOuterSubtileMap(4);
        outerMaxX = outerMinX + outerWidthX - 1;
        outerMaxY = outerMinY + outerHeightY - 1;
        %define the overlap (boundary region)
        overlapX = innerMinX - outerMinX;
        overlapY = innerMinY - outerMinY;
        %define the offset from the local subtile to the global coords of the image
        localWindoxOffsetX = outerMinX - 1;
        localWindoxOffsetY = outerMinY - 1;

        % get the subimage
        Isub = I(outerMinY:outerMaxY,outerMinX:outerMaxX); %crop real image to the current subimage only for fitting
        maxIntensity = max(Isub(:)); %max intensity in the subtile
        bgEst = mean(Isub(:)); %initial guess background

        %get data points with cluster centers within the entire
        %subimage (note that clusters in the border region are only included for improved fit--
        %only the values from clusters within the inner window will be reported)
        outerIdx = find(roundedAdjY >= outerMinY & roundedAdjX >= outerMinX & roundedAdjY <= outerMaxY & roundedAdjX <= outerMaxX);
        numInOuter = length(outerIdx); %number of points within the subimage

        % initial guess x/y position is taken from the registered cluster
        % center; adjust from global image to local subimage coordinates
        yEst = adjY(outerIdx) - localWindoxOffsetY;
        xEst = adjX(outerIdx) - localWindoxOffsetX;
        
        %get indices of the fits that fall within the inner fit window.
        localInnerWindowMinX = overlapX;
        localInnerWindowMinY = overlapY;
        localInnerWindowMaxX = outerWidthX - overlapX + 1;
        localInnerWindowMaxY = outerHeightY - overlapY + 1;
        innerIdx = find(xEst > localInnerWindowMinX & yEst > localInnerWindowMinY & xEst <= localInnerWindowMaxX & yEst <= localInnerWindowMaxY); %get indices of only the fits that fall within the inner fit window
        numInInner = length(innerIdx); %number of points within the subimage

        % determine if any cluster in the frame overlaps with junk regions
        % cleaned by CleanImage
        cPos = sub2ind(size(I),roundedAdjY(outerIdx),roundedAdjX(outerIdx));
        idx2 = ismember(cPos,zeroPix);
        uncleaned = sum(idx2)==0; %true if none of the regions in this subtile have been cleaned by CleanImage

        % if clusters exist in this region and none of them overlap with junk regions...
        if numInInner>0 && uncleaned 
            % initialize fit parameter arrays
            nsz = 4*numInOuter+1; %fit parameter arrays must be variable length because each iteration may be fitting a different number of clusters
                                  %each cluster will have 4 parameters: amplitude, sigma, x position, & y position
                                  %there is also one global background param

            c0 = zeros(1,nsz); %initial guess array
            ub = zeros(1,nsz); %upper bound array
            lb = zeros(1,nsz); %lower bound array

            % initial guess amplitude
            aEst = Isub(sub2ind(size(Isub),round(yEst),round(xEst))); %use the closest pixel value to the cluster center as a first guess
            aEst = max((aEst-bgEst),1); %adjust estimate for background

            % parameters and bounds, reassigned for each tile due to variable numbers of tiles
            ampPos = 1:4:nsz-1;
            xCentPos = 2:4:nsz-1;
            yCentPos = 3:4:nsz-1;
            sigPos = 4:4:nsz-1;

            % amplitudes
            lb(ampPos) = 0;  
            c0(ampPos) = aEst;
            ub(ampPos) = maxIntensity*1.2;
            
            % x center positions
            lb(xCentPos) = xEst - GlobalVars.positionTolerance;
            c0(xCentPos) = xEst;
            ub(xCentPos) = xEst + GlobalVars.positionTolerance;
            
            % y center positions
            lb(yCentPos) = yEst - GlobalVars.positionTolerance;
            c0(yCentPos) = yEst;
            ub(yCentPos) = yEst + GlobalVars.positionTolerance;
            
            % sigmas
            lb(sigPos) = GlobalVars.multiGauss_sigmaRange(1);  
            c0(sigPos) = GlobalVars.multiGauss_sigmaRange(2);
            ub(sigPos) = GlobalVars.multiGauss_sigmaRange(3);
            
            % global background
            lb(end) = 0;
            c0(end) = bgEst;
            ub(end) = maxIntensity;
             
            % set up the local y,x coordinate matrix for fitting
            currFitWindowY = fitWindowY(1:outerWidthX,1:outerHeightY);
            currFitWindowX = fitWindowX(1:outerWidthX,1:outerHeightY);
            yx = zeros(outerHeightY*outerWidthX,2);
            yx(:,1) = currFitWindowY(:);
            yx(:,2) = currFitWindowX(:);

            % fit n 2D gaussians to the subimage with least squares
            transposeIsub = Isub';
            cc = lsqcurvefit(fun,c0,yx,transposeIsub(:),lb,ub,options);

            %### DEBUG ###
            if(GlobalVars.debug)
                %save boundary region data fits for debugging visualization
                cDat.fitAmplitudes(unexcisedIndices(outerIdx)) = cc(ampPos);  % amplitudes
                cDat.fitY(unexcisedIndices(outerIdx)) = cc(yCentPos) + localWindoxOffsetY;  % y center positions
                cDat.fitX(unexcisedIndices(outerIdx)) = cc(xCentPos) + localWindoxOffsetX;  % x center positions
                cDat.fitSigmas(unexcisedIndices(outerIdx)) = cc(sigPos);  % sigmas
                cDat.alreadyFit(unexcisedIndices(outerIdx)) = false;
            end
            %### END DEBUG ###
            
            % output fit parameters for clusters within the inner fit
            % window to the ClusterData, indexed appropriately
            cDat.fitAmplitudes(unexcisedIndices(outerIdx(innerIdx))) = cc(ampPos(innerIdx));  % amplitudes
            cDat.fitY(unexcisedIndices(outerIdx(innerIdx))) = cc(yCentPos(innerIdx)) + localWindoxOffsetY;  % y center positions
            cDat.fitX(unexcisedIndices(outerIdx(innerIdx))) = cc(xCentPos(innerIdx)) + localWindoxOffsetX;  % x center positions
            cDat.fitSigmas(unexcisedIndices(outerIdx(innerIdx))) = cc(sigPos(innerIdx));  % sigmas
            cDat.alreadyFit(unexcisedIndices(outerIdx(innerIdx))) = true;
            bgImage(innerMinY:innerMaxY,innerMinX:innerMaxX) = cc(end); %update background image

            %### DEBUG ###
            %plot fits individual subtile by subtile for debugging
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:QuantifyImage: Would you like to see the multi-gaussian fit for the full current subtile?');
                if(showDebug)
                    if(~exist('figureHandle','var') || ~ishandle(figureHandle))
                        figureHandle = figure();
                    end
                    figureHandle = PlotFun.plotThreePanelFit(I, cDat, cc(end), currOuterSubtileMap, currOuterSubtileMap, [], [], figureHandle);
                end
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:QuantifyImage: Would you like to see the multi-gaussian fit for the inner window only?');
                if(showDebug)
                    if(~exist('figureHandle','var') || ~ishandle(figureHandle))
                        figureHandle = figure();
                    end
                    figureHandle = PlotFun.plotThreePanelFit(I, cDat, cc(end), currOuterSubtileMap, currInnerSubtileMap, [], GlobalVars.overlap, figureHandle);
                end
            end
            %### END DEBUG ###
        else
            bgImage(innerMinY:innerMaxY,innerMinX:innerMaxX) = bgEst; %if this subtile was not fit, just set the background to the background estimate for the tile
        end
        
        %progress bar
        progbar(subtileIdx,1);
    end
    
    %### DEBUG ###
    if(GlobalVars.debug)
        showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:QuantifyImage: Would you like to see the execution time profile for the fitting routine?');
        if(showDebug)    
            profile viewer;
            profile off;
        end
    end
    % ### END DEBUG ###
end
