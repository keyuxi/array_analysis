% Class for registering sequencing data to cluster images

classdef RegistrationFun
     methods (Static)
         
        % Global registration of a central tile of the image with all seq data
        function [cDat, globalOffsetY, globalOffsetX, success] = RegisterGlobal(image, cDat)
            %define the central tile to register (size set in GlobalVars)
            centralSubtileCoords = ImageFun.getCentralSubtileCoords(size(image));

            %define the bounds of the data to register (the entire set of data)
            dataBounds = [min(cDat.currX()) min(cDat.currY()) max(cDat.currX()) max(cDat.currY())];
            
            %register
            cDat.resetAlreadyMoved();
            %find if GlobalVars has offset limits for global registration
            if(sum(ismember(properties(GlobalVars),'maxOffsetGlobalRegistration'))==1)
                [cDat, globalOffsetY, globalOffsetX, success, correlationPeakHeight] = RegistrationFun.RegisterTile(cDat, image, centralSubtileCoords, dataBounds, GlobalVars.maxOffsetGlobalRegistration);
            else
                [cDat, globalOffsetY, globalOffsetX, success, correlationPeakHeight] = RegistrationFun.RegisterTile(cDat, image, centralSubtileCoords, dataBounds);
            end
            
            if(~success) %if the global registration failed
                error('RegisterGlobal:RegisterGlobal','Global registration failed (peak height = %f, threshold for correlation success = %f).  Analysis cannot continue. Is this the correct tile? Was the correct data source/scaling factor chosen?',correlationPeakHeight,GlobalVars.correlationSuccessPeakHeight);
            end
            
            cDat.updateExcised(~cDat.alreadyMoved); %excise anything that didn't get moved during global registation
        end % END RegisterGlobal

        
        
        % Break up the camera image into subtiles, and register each subtile
        % individually. This helps to handle spherical abberations. The maximum
        % allowed x and y offsets are passed in, so that this procedure is
        % constrained to not have each subtile registration stray too far from the
        % globally-registered values. (The clusters_x and clusters_y data points
        % are assumed to have been adjusted to account for the global registration
        % already.)
        function [cDat, numSuccessful] = RegisterImageInSubtiles(cDat, rawImage, subtileMap, maxOffsetDeltaY, maxOffsetDeltaX)
            numSuccessful = 0;
            dataBounds = zeros(1,4);
            
            cDat.resetAlreadyMoved();
            for i = 1:size(subtileMap,1)
                subImageRect = subtileMap(i,:);

                % define the bounds of the data over this subtile, based on the current cluster positions (which should at least have been globally registered already)
                dataBounds(1) = subImageRect(1);
                dataBounds(2) = subImageRect(2);
                dataBounds(3) = dataBounds(1) + subImageRect(3);
                dataBounds(4) = dataBounds(2) + subImageRect(4);
                if ~isempty(find(cDat.adjY >= dataBounds(2) & cDat.adjY < dataBounds(4) & cDat.adjX >= dataBounds(1) & cDat.adjX < dataBounds(3) & ~cDat.alreadyMoved, 1))
                    % if one or more clusters were located in this subtile, then proceed with subtile registration run registration on this subtile, with only the cluster data 
                    % that lies on top of the subtile pre-registration (within dataBounds)  This makes the assumption that global registration has brought us close)
                    [cDat, localOffsetY, localOffsetX, currTileRegistrationSuccess] = RegistrationFun.RegisterTile(cDat, rawImage, subImageRect, dataBounds, [maxOffsetDeltaY maxOffsetDeltaX]);
                    if(currTileRegistrationSuccess)
                        numSuccessful = numSuccessful + 1;
                    end
                end
            end
            cDat.updateExcised(~cDat.alreadyMoved); %excise anything that didn't get moved during this round of registration
        end % END RegisterImageInSubtiles
        
        
        
        % Register a part of a camera image (a "subtile") to a set of cluster
        % positions. Performs the registration subject to the constraint that the
        % offset between the camera subtile image and the cluster data cannot be
        % larger than the maxOffset specified. 
        function [cDat, offsetY, offsetX, registrationSuccess, correlationPeakHeight] = RegisterTile(cDat, rawImage, subImageRect, dataBounds, maxOffset)
            %### DEBUG ###
            if(GlobalVars.debug)
                global debugFigMan;
            end
            %### END DEBUG ###
            
            if(~exist('maxOffset', 'var'))
                maxOffset = [Inf Inf]; %register with a maxiumum x,y offset from the center of the correlation matrix ([Inf Inf] is an unbounded search)
            elseif((maxOffset(1) < 1) || (maxOffset(2) < 1))
                error('RegisterTile:RegisterTile','max offset for bounded correlation must be at least 1x1 pixel');
            end

            % crop out a subtile region of the real image for correlation
            subtileImage = imcrop(rawImage, subImageRect);
            subtileOffset = [(1 - subImageRect(2)) (1 - subImageRect(1))]; %the offset of the subimage relative to the raw image

            % Generate synthetic image from cluster data for cross-correlation against the real image

            % find all the points that fall within the specified data bounds
            subtileIndices = (cDat.adjY >= dataBounds(2) & cDat.adjY < dataBounds(4) & cDat.adjX >= dataBounds(1) & cDat.adjX < dataBounds(3) & ~cDat.alreadyMoved & ~cDat.excised);
            
            if(size(find(subtileIndices & cDat.currSubset),1) < GlobalVars.minNumClustersForRegistration) %if no or too few valid clusters were found in this subtile..
                %do not adjust cluster positions if registration failed
                registrationSuccess = false;
                correlationPeakHeight = 0;
                offsetY = 0;
                offsetX = 0;
                %excise any other clusters of the current subtile, because registration has failed
                cDat.updateExcised(subtileIndices);
                return
            end
            
            %only use points in the subtile and that are part of the
            %current subset (as defined by filter matching) for registration
            subtileAdjY = cDat.adjY(subtileIndices & cDat.currSubset);
            subtileAdjX = cDat.adjX(subtileIndices & cDat.currSubset);
            
            % define the bounds of the synthetic image such that they are guaranteed
            % to be larger than the central region of the real image (normxcorr2 demands that
            % the first argument (real image) have smaller dimensions than the second argument(synth image))
            synthImageSizeY = max( ceil(max(subtileAdjY)-min(subtileAdjY))+1 , subImageRect(4)+2 );
            synthImageSizeX = max( ceil(max(subtileAdjX)-min(subtileAdjX))+1 , subImageRect(3)+2 );

            %create a synthetic image for cross-correlation
            synthImg = ImageFun.makeSynthImagePoints(subtileAdjY, subtileAdjX, synthImageSizeY, synthImageSizeX, GlobalVars.synthIntensityValue);

            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG: AnalyseImage: Would you like to see the synthetic image of the sequence data positions?');
                if(showDebug)
                    PlotFun.plotClusterImage(synthImg, 'synthetic image of sequencing data for cross-correlation');
                end
                drawnow;
            end
            %### END DEBUG ###

            % Normalized crosscorrelation of the subtile image versus the synthetic image
            corrMatrix = normxcorr2(subtileImage, synthImg);

            %fit the correlation matrix to a gaussian to find the peak
            %correlation to sub-pixel resolution.  The correlation peak height is
            %used as a measure of registration success
            [registrationSuccess yPeak xPeak correlationPeakHeight] = RegistrationFun.fitCorrelationGaussian(corrMatrix, maxOffset);

            if(registrationSuccess) %if the cross-correlation for this tile was successful
                corrOffset = [(yPeak - size(subtileImage,1)) (xPeak - size(subtileImage,2))]; %the offset of the correlation peak relative to the subtile image size

                % sum the relative offsets of the position of the central subimage (synthetic data) registered onto the actual image
                offsetY = min(subtileAdjY) + corrOffset(1) + subtileOffset(1) - 1;
                offsetX = min(subtileAdjX) + corrOffset(2) + subtileOffset(2) - 1;

                % adjust cluster positions by offset
                cDat.adjY(subtileIndices) = cDat.adjY(subtileIndices) - offsetY;
                cDat.adjX(subtileIndices) = cDat.adjX(subtileIndices) - offsetX;
                cDat.alreadyMoved(subtileIndices) = true;
            else
                %do not adjust cluster positions if registration failed
                offsetY = 0;
                offsetX = 0;
                %excise the clusters of the currernt subtile that on which registration has failed
                cDat.updateExcised(subtileIndices);
            end
        end % END RegisterTile

        
        
        %fit the peak of a crosscorrelation matrix to a 2D gaussian to get
        %sub-pixel resolution
        function [registrationSuccess yPeak xPeak correlationPeakHeight] = fitCorrelationGaussian(corrMatrix, maxOffset)
            %### DEBUG ###
            if(GlobalVars.debug)
                global debugFigMan;
            end
            %### END DEBUG ###
            
            if(~exist('maxOffset','var'))
                maxOffset = [Inf Inf]; %maxiumum x,y offset from the center of the correlation matrix ([Inf Inf] is an unbounded search)
            end

            registrationSuccess = false;

            %size of the correlation matrix
            sizeCorrMatrix = size(corrMatrix);

            %boolean variable - the search is considered bounded if any non-infinite bound is passed in
            bounded = ((maxOffset(1) ~= Inf) && (maxOffset(2) ~= Inf)); 

            %...calculate the bounds in which we will look for a correlation peak
            if(bounded) %if the search is bounded confine the search to the bounded distance from the center of the matrix
                centerPositionY = round(sizeCorrMatrix(1)/2);
                centerPositionX = round(sizeCorrMatrix(2)/2);
                minY = max(1,centerPositionY-maxOffset(1));
                maxY = min(sizeCorrMatrix(1),centerPositionY+maxOffset(1));
                minX = max(1,centerPositionX-maxOffset(2));
                maxX = min(sizeCorrMatrix(2),centerPositionX+maxOffset(2));
            else %if unbounded set global bounds
                minY = 1;
                maxY = sizeCorrMatrix(1);
                minX = 1;
                maxX = sizeCorrMatrix(2);
            end

            putativePeak = zeros(1,2);
            %...assign the putative peak based on a peak finder
            if(bounded)%if bounded...
                %...use the max within the bounded window
                boundedCorrMatrix = corrMatrix(minY:maxY,minX:maxX); %get the bounded correlation matrix
                [putativePeak(1) putativePeak(2) correlationPeakHeight] = RegistrationFun.findCorrelationPeak(boundedCorrMatrix);
                putativePeak(1) = putativePeak(1) + minY - 1; %adjust for the offset of the bounded window from the global correlation matrix
                putativePeak(2) = putativePeak(2) + minX - 1;
            else %if unbounded...
                %...use the global correlation matrix
                [putativePeak(1) putativePeak(2) correlationPeakHeight] = RegistrationFun.findCorrelationPeak(corrMatrix);
            end

            %get a small window around the putative peak to fit to a gaussian
            fitWindowYmin = max(putativePeak(1)-GlobalVars.corrPeak_window,minY);
            fitWindowYmax = min(putativePeak(1)+GlobalVars.corrPeak_window,maxY);
            fitWindowXmin = max(putativePeak(2)-GlobalVars.corrPeak_window,minX);
            fitWindowXmax = min(putativePeak(2)+GlobalVars.corrPeak_window,maxX);
            
            %define the data window to be fit
            corrPeakWindow = corrMatrix(fitWindowYmin:fitWindowYmax,fitWindowXmin:fitWindowXmax);
            sizeCorrPeakWindow = size(corrPeakWindow);
            [corrPeakWindowX corrPeakWindowY] = meshgrid(1:sizeCorrPeakWindow(2),1:sizeCorrPeakWindow(1));
            yx(:,1) = corrPeakWindowY(:) + fitWindowYmin - 1;
            yx(:,2) = corrPeakWindowX(:) + fitWindowXmin - 1;
                
            %                                   A                                      y0                                             x0                      sigma                        bg
            lb           = [                    0, putativePeak(1)-GlobalVars.corrPeak_peakBound, putativePeak(2)-GlobalVars.corrPeak_peakBound, GlobalVars.corrPeak_sigmaRange(1), min(corrMatrix(:))];
            initialGuess = [correlationPeakHeight,                               putativePeak(1),                               putativePeak(2), GlobalVars.corrPeak_sigmaRange(2), min(corrMatrix(:))];
            ub           = [                    1, putativePeak(1)+GlobalVars.corrPeak_peakBound, putativePeak(2)+GlobalVars.corrPeak_peakBound, GlobalVars.corrPeak_sigmaRange(3), max(corrMatrix(:))];
            
            options = optimset( 'TolX', GlobalVars.corrPeak_TolX,... % TolX = convergence criteron on the change in parameters one iteration to the next
                                'TolFun', GlobalVars.corrPeak_TolFun,... % TolX = convergence criteron on the change in the value of the objective function one iteration to the next
                                'MaxFunEvals', GlobalVars.corrPeak_MaxFunEvals,... %max number of times objective function can be evaluated
                                'Display', 'off',... %display no output from optimization
                                'jacobian', 'on',... %use the jacobian
                                'DerivativeCheck', 'off',... %do not compare user-supplied jacobian (gradient of objective) to calculated finite-differencing derivatives
                                'LargeScale', 'on');

            fun = @FitFun.gaussian_2D;

            [fitParams, resnorm] = lsqcurvefit(fun,initialGuess,yx,corrPeakWindow(:),lb,ub,options);
            
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:fitCorrelationGaussian: Would you like to see the correlation matrix 2D gaussian fit?');
                if(showDebug)
                    msg = sprintf('DEBUG:RegistrationFun:fitCorrelationGaussian: fit correlation peak.  Peak height = %f', correlationPeakHeight);
                    PlotFun.plotFitCorrelationMatrix(corrPeakWindow,fitWindowYmin,fitWindowXmin,fitParams,msg);
                end
                drawnow;
            end
            %### END DEBUG ###
            
            if(correlationPeakHeight > GlobalVars.correlationSuccessPeakHeight) 
                registrationSuccess = true;
            end

            yPeak = fitParams(2);
            xPeak = fitParams(3);
        end % END fitCorrelationGaussian
        
        

        %find the peak of a crosscorrelation matrix
        function [yPeak xPeak correlationPeakHeight] = findCorrelationPeak(corrMatrix, maxPeakWidth)
            %### DEBUG ###
            if(GlobalVars.debug)
                global debugFigMan;
            end
            %### END DEBUG ###
            
            % find if maxPeakWidth is set
            if(~exist('maxPeakWidth','var'))
                % if maxPeakWidth is set in globalvars, use that value. Ohterwise set it to ten.
                if(sum(ismember(properties(GlobalVars),'corrPeak_background_window'))==1);
                    maxPeakWidth = GlobalVars.corrPeak_background_window;
                else
                    maxPeakWidth = 10; %m ax expected width of the correlation peak, in pixels
                end
            end

            sizeCorrMatrix = size(corrMatrix);
            bgCorrMatrix = ImageFun.getBackgroundFromClusterImage(corrMatrix, maxPeakWidth);
            bgSubtractedCorrMatrix = corrMatrix - bgCorrMatrix;       
    
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:findCorrelationPeak: Would you like to see the correlation matrix subtraction for peak finding?');
                if(showDebug)
                    PlotFun.plotCorrelationMatrix(corrMatrix,1,1,'before bg subtraction');
                    PlotFun.plotCorrelationMatrix(bgCorrMatrix,1,1,'background matrix');
                    PlotFun.plotCorrelationMatrix(bgSubtractedCorrMatrix,1,1,'after bg subtraction');
                end
                drawnow;
            end
            %### END DEBUG ###

            [currPeakMax, imax] = max(bgSubtractedCorrMatrix(:)); %find the max
            [yPeak xPeak] = ind2sub(sizeCorrMatrix, imax);

            boundary = floor(maxPeakWidth/2);

            currBoxMax = 0;
            boxPadding = 2;
            
            YnegBound = max((yPeak-boundary), 1);
            YposBound = min((yPeak+boundary), sizeCorrMatrix(1));
            XnegBound = max((xPeak-boundary), 1);
            XposBound = min((xPeak+boundary), sizeCorrMatrix(2));
            
            if((yPeak-boundary) > boxPadding)
                currBoxMax = max([currBoxMax bgSubtractedCorrMatrix((yPeak-boundary),XnegBound:XposBound)]);
            end
            if((yPeak+boundary) <= (sizeCorrMatrix(1)-boxPadding))
                currBoxMax = max([currBoxMax bgSubtractedCorrMatrix((yPeak+boundary),XnegBound:XposBound)]);
            end
            if((xPeak-boundary) > boxPadding)
                currBoxMax = max([currBoxMax bgSubtractedCorrMatrix(YnegBound:YposBound,(xPeak-boundary))']);
            end
            if((xPeak+boundary) <= (sizeCorrMatrix(2)-boxPadding))
                currBoxMax = max([currBoxMax bgSubtractedCorrMatrix(YnegBound:YposBound,(xPeak+boundary))']);
            end

            correlationPeakHeight = (currPeakMax-currBoxMax);
        end % END findCorrelationPeak

        
        
        %apply a pre-computed registration offset map to a set of clusters
        function [cDat,localOffsets_y,localOffsets_x] = ApplyRegistrationOffset(rawImage,cDat,currRegistrationOffsetMap,imageSize,maskImage,globalOffset_y,globalOffset_x)
            disp('Applying registration offset map...');
            %If the global registration offset is passed in, it implies the global
            %offset has NOT ALREADY BEEN APPLIED to the adjusted cluster
            %positions
            if(exist('globalOffset_y','var'))
                cDat.adjY = cDat.adjY - globalOffset_y; %apply the global offset
            end
            if(exist('globalOffset_x','var'))
                cDat.adjX = cDat.adjX - globalOffset_x; %apply the global offset
            end
            
            offsetMapParams_y = currRegistrationOffsetMap(1,:);
            offsetMapParams_x = currRegistrationOffsetMap(2,:);

            %calculate the local offsets from the offset map
            localOffsets_y = FitFun.quadraticSurface(offsetMapParams_y(1:end-2),[cDat.adjY cDat.adjX],offsetMapParams_y(end-1),offsetMapParams_y(end));
            localOffsets_x = FitFun.quadraticSurface(offsetMapParams_x(1:end-2),[cDat.adjY cDat.adjX],offsetMapParams_y(end-1),offsetMapParams_y(end));

            %apply local offsets
            cDat.adjY = cDat.adjY + localOffsets_y;
            cDat.adjX = cDat.adjX + localOffsets_x;

            %excise the data to remove sequences that are out of bounds
            [cDat] = RegistrationFun.ExciseSeqData(cDat,imageSize,maskImage);

        end % END ApplyRegistrationOffset

        
        
        % Hierarchical registration of cluster data to camera image
        % Strategy:
        %   1. Perform a "global" registration, i.e. a registration of a large,
        %   central chunk of the camera image against all the cluster data. Adjust
        %   the cluster data using this offset.
        %
        %   2. Divide the image into "subtiles", e.g. 3x3=9 different subimages,
        %   and register each subtile separately against the cluster points that
        %   lie within the subtile based on the global registration. Adjust the
        %   cluster values for each subtile based on the new offset found for that
        %   subtile. The registration offsets are constrained to be in some range
        %   (e.g. +/- 10 pixels); this prevents the subtile registration from
        %   finding a new and completely wrong registration, which might happen if
        %   the number of clusters within a single subtile is small.
        %
        %   3. Divide the image into even smaller subtiles (e.g. 4x4=16), and
        %   perform registration on subtiles, as in step 2, except now the starting
        %   point is the set of adjusted cluster positions at the output of step 2.
        %   Typically the constraints on how many pixels the offset can be here
        %   will be smaller than in stage 2.
        %
        %   4. Repeat step 3 until desired level of subtile size is reached to
        %   handle the degree of image distortion being dealt with. (For GA imaging
        %   station, versus MiSeq cluster data, it seems that just two levels of
        %   subtile registration is sufficient.)
        function [cDat, localOffsetsY, localOffsetsX, globalOffsetY, globalOffsetX] = RegisterHierarchical(image, cDat, filterColormap)
            %### DEBUG ###
            if(GlobalVars.debug)
                global debugFigMan;

                if~exist('filterColormap','var')
		    filterColormap = containers.Map();
		end
            end
            %### END DEBUG ###
            
            disp('Hierarchical registration:');
            % First perform global registration of a central tile of the image with all seq data
            disp('    ...initial global registration...');
            [cDat, globalOffsetY, globalOffsetX, success] = RegistrationFun.RegisterGlobal(image, cDat);
            if(~success) %if the global registration failed
                error('RegisterHierarchical:RegisterHierarchical','Global registration failed.  Analysis cannot continue. Is this the correct tile? Was the correct data source/scaling factor chosen?  The Global variable "correlationSuccessPeakHeight" can be changed to adjust stringency, if desired.');
            end
            
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:RegisterHierarchical: Would you like to see the initial global registration?');
                if(showDebug)
                      PlotFun.plotGlobalRegistrationByFilter(image, cDat, [], [], 'initial global registration', filterColormap);
                end
                drawnow;
                debugFigMan.clearLocal(); %clear debug dismissals from the previous level of registration
            end
            %### END DEBUG ###

            % now register subtiles of image against all seq data, subject to
            % constraint of the offsets we obtained from the global image
            % registration
            nTiles = GlobalVars.numTilesCoarseHierarchicalRegistration;
            totalNumTiles = nTiles(1)*nTiles(2);
            maxOffset = GlobalVars.maxOffsetCoarseHierarchicalRegistration;
            disp(['    ...coarse hierarchical subtile registration (' num2str(nTiles(1)) 'x' num2str(nTiles(2)) ')...']);
            subtileMap = ImageFun.getSubtileMapCoords(size(image), nTiles(1), nTiles(2));
            [cDat, numSuccessful] = RegistrationFun.RegisterImageInSubtiles(cDat, image, subtileMap, maxOffset(1), maxOffset(2));
            disp(['            ' num2str(numSuccessful) '/' num2str(totalNumTiles) ' subtiles successfully registered']);
            if(numSuccessful <= 0) %if the coarse hierarchical registration failed
                error('RegisterHierarchical:RegisterHierarchical','Coarse hierarchical subtile registration failed.  Analysis cannot continue. Is this the correct tile? Was the correct data source/scaling factor chosen?');
            end
            
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:RegisterHierarchical: Would you like to see the coarse subtile registration?');
                if(showDebug)
                    PlotFun.plotGlobalRegistrationByFilter(image, cDat, [], [], 'coarse subtile registration', filterColormap);
                end
                drawnow;
                debugFigMan.clearLocal(); %clear debug dismissals from the previous level of registration
            end
            %### END DEBUG ###
            
            % register again using smaller subtiles
            nTiles = GlobalVars.numTilesFineHierarchicalRegistration;
            totalNumTiles = nTiles(1)*nTiles(2);
            maxOffset = GlobalVars.maxOffsetFineHierarchicalRegistration;
            disp(['    ...fine hierarchical subtile registration (' num2str(nTiles(1)) 'x' num2str(nTiles(2)) ')...']);
            subtileMap = ImageFun.getSubtileMapCoords(size(image), nTiles(1), nTiles(2));
            [cDat, numSuccessful] = RegistrationFun.RegisterImageInSubtiles(cDat, image, subtileMap, maxOffset(1), maxOffset(2));
            disp(['            ' num2str(numSuccessful) '/' num2str(totalNumTiles) ' subtiles successfully registered']);
            if(numSuccessful <= 0) %if the fine hierarchical registration failed
                error('RegisterHierarchical:RegisterHierarchical','Fine hierarchical subtile registration failed.  Analysis cannot continue. Is this the correct tile? Was the correct data source/scaling factor chosen?');
            end
            
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:RegisterHierarchical: Would you like to see the fine subtile registration?');
                if(showDebug)
                    PlotFun.plotGlobalRegistrationByFilter(image, cDat, [], [], 'fine subtile registration', filterColormap);
                end
                drawnow;
                debugFigMan.clearLocal();
            end
            %### END DEBUG ###
            
            % calculate "local" offsets (i.e. offsets used per-point, relative to
            % the global offset calculated from the first registration)
            localOffsetsY = globalOffsetY + cDat.adjY - cDat.origY;
            localOffsetsX = globalOffsetX + cDat.adjX - cDat.origX;
        end % END RegisterHierarchical
         

        
        % Remove all cluster points that do not lie within the "valid" part of
        % the image. The GA imaging station image is of a smaller region than
        % the MiSeq images, so already there are many cluster points that fall
        % outside the image, and these should be excised. Also remove any points 
        % that were marked for deletion because their subtile could not register.
        % Tiles are marked for deletion by setting their position to NaN,NaN

        % A mask image may also be used to explicitly remove data that align to
        % a particular region--The mask image is defined so that pixels with value
        % 255 are defined as being within the usable field-of-view, and pixels with
        % value 0 are not. All sequencing data points that fall on top of mask
        % image pixels with value 255 are returned, and all others are dropped.

        % Input cluster positions should have been scaled  to be in units of image pixels, and should have been registered to the camera image.
        function [cDat] = ExciseSeqData(cDat, imageSize, maskImage)
            % Exclude all clusters that fall outside of the image
            % boundaries (also exclude an additional small boundary around
            % the image that cannot be fit well due to edge effects)

            out_of_bounds = (round(cDat.adjX) <= GlobalVars.overlap | round(cDat.adjX) >= (imageSize(2)-GlobalVars.overlap) | round(cDat.adjY) <= GlobalVars.overlap | round(cDat.adjY) >= (imageSize(1)-GlobalVars.overlap));
            cDat.updateExcised(out_of_bounds);

            % exclude all clusters that fall within disallowed regions of the mask
            if(~isempty(maskImage))
                maskImageSize = size(maskImage);
                if(~isequal(maskImageSize,imageSize))
                    errMsg = ['Mask image (' num2str(maskImageSize(1)) 'x' num2str(maskImageSize(2)) ') must be the same size as the image being analyzed (' num2str(imageSize(1)) 'x' num2str(imageSize(2)) ')'];
                    error('RegistrationFun:ExciseSeqData',errMsg);
                end
                                
                % Exclude all clusters that don't fall within allowed part of mask
                unexcisedIndices = find(~cDat.excised);
                unexcisedAdjY = cDat.adjY(unexcisedIndices);
                unexcisedAdjX = cDat.adjX(unexcisedIndices);
                out_of_mask = (maskImage(sub2ind(size(maskImage),round(unexcisedAdjY),round(unexcisedAdjX))) ~= 255); % for each x,y cluster position, find the closest pixel in the mask image, and if the mask pixel value is 255, include the index of that x,y cluster in the list of valid points
                toBeExcised = false(size(cDat.excised));
                toBeExcised(unexcisedIndices(out_of_mask)) = true;
                cDat.updateExcised(toBeExcised);
            end
        end
        
        
        %calcuate the registration offset map from an image with
        %hierarchical registration
        function [registrationOffsetMap, globalOffset_y, globalOffset_x] = calculateRegistrationOffsetMap(image, cDat, maskImage, filterColormap)
            %### DEBUG ###
            if(GlobalVars.debug)
                global debugFigMan;

                if~exist('filterColormap','var')
                    filterColormap = containers.Map();
                end
            end
            %### END DEBUG ###

            % remove background from the image with morphological opening
            % (this greatly improves cross correlation performance)
            bgImage = ImageFun.getBackgroundFromClusterImage(image);
            bgSubtractedImage = image - bgImage;

            %hierarchical registration
            disp('Computing the registration offset map with hierarchical registration...')
            [cDat, localOffsets_y, localOffsets_x, globalOffset_y, globalOffset_x] = RegistrationFun.RegisterHierarchical(bgSubtractedImage, cDat, filterColormap);

            % Remove all cluster points that do not lie within the "valid" part of
            % the image. The GA imaging station image is of a smaller region than
            % the MiSeq images, so already there are many cluster points that fall
            % outside the image, and these should be excised. Also remove any points 
            % that were marked for deletion because their subtile could not register
            % a mask image may also be used to explicitly remove data that align to
            % a particular region
            [cDat] = RegistrationFun.ExciseSeqData(cDat, size(image), maskImage);
            
            %get a position map of registration offsets relative to the global registration
            disp('Fitting the map to the hierarchicaly registered clusters...')
            offsetMapParams_y = PlotFun.fitRegistrationOffset(cDat.currY(), cDat.currX(), localOffsets_y(cDat.currSubset));
            offsetMapParams_x = PlotFun.fitRegistrationOffset(cDat.currY(), cDat.currX(), localOffsets_x(cDat.currSubset));

            registrationOffsetMap = [offsetMapParams_y; offsetMapParams_x];
                        
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:calculateRegistrationOffsetMap: Would you like to see final registered data, post excision?');
                if(showDebug)
%                     PlotFun.plotRegistration(image, [cDat.currY() cDat.currX()], min(image(:)), max(image(:))-min(image(:)), 'registered data, post-excision');
                    PlotFun.plotGlobalRegistrationByFilter(image, cDat, min(image(:)), max(image(:))/GlobalVars.goldenMean, 'registered data, post-excision', filterColormap);
                end
                drawnow;
            end

            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:calculateRegistrationOffsetMap: Would you like to see registration offset map fits?');
                if(showDebug)
                    PlotFun.plotRegistrationOffsetFit(cDat.currY(),cDat.currX(),localOffsets_y(cDat.currSubset),'y registration offsets, fit',offsetMapParams_y);
                    PlotFun.plotRegistrationOffsetFit(cDat.currY(),cDat.currX(),localOffsets_x(cDat.currSubset),'x registration offsets, fit',offsetMapParams_x);
                end
                drawnow;
            end
            %### END DEBUG ###
            
        end % END calculateRegistrationOffsetMap
        
        
        
        function saveRegistrationOffsetMap(registrationOffsetMap, filename)
            save(filename,'registrationOffsetMap','-ascii');
        end
        
        
        
        function registrationOffsetMap = loadRegistrationOffsetMap(filename)
            registrationOffsetMap = load(filename);
        end
        
        
        %generate a field of view mask, in particular for images with a
        %ciruclar field of view on a square image that have unilluminated
        %regions in the corner.  Creates a mask image that will be overlaid
        %on top of the data. Regions with value 0 will be discarded and
        %regions with value 255 (max) will be kept.
        function maskImage = autoGenerateFOVmask(image, maskFilename)
            %### DEBUG ###
            if(GlobalVars.debug)
                global debugFigMan;
            end
            %### END DEBUG ###
            
            maskImage = zeros(size(image)); %initialize mask image array

            %get background illumination
            bgImage = ImageFun.getBackgroundFromClusterImage(image); 
            
            minVal = min(bgImage(:));
            
            %subtract minimum value
            bgImage = bgImage - minVal;
            
            %smooth background image
            smoothedImage = imfilter(bgImage,fspecial('disk', GlobalVars.smoothingFilterRadius),'replicate');           
            
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:autoGenerateFOVmask: Would you like to see debugging figures for FOV mask generation?');
                if(showDebug)
                    PlotFun.plotClusterImage(image, 'DEBUG:RegistrationFun:autoGenerateFOVmask: Starting with the Raw Image');
                    PlotFun.plotImageIntensitySurface(image, 'DEBUG:RegistrationFun:autoGenerateFOVmask: Starting with the Raw Image');
                    PlotFun.plotClusterImage(smoothedImage, 'DEBUG:RegistrationFun:autoGenerateFOVmask: Smoothed background illumination');
                    PlotFun.plotImageIntensitySurface(smoothedImage, 'DEBUG:RegistrationFun:autoGenerateFOVmask: Smoothed background illumination');
                end
                drawnow;
            end
            %### END DEBUG ###
            
            %determine the cutoff between illuminated and unilluminated
            %regions as a fraction of the mean intensity
            meanIntensity = mean(smoothedImage(:));
            illuminatedRegionCutoff = meanIntensity*GlobalVars.maskIntensityCutoff;
            
            %find illuminated image regions
            illuminatedIndices = find(smoothedImage>=illuminatedRegionCutoff);
            
            %set illuminated image regions to the max value            
            maskImage(illuminatedIndices) = 255;
            
            %### DEBUG ###
            if(GlobalVars.debug)
                showDebug = debugFigMan.debugStop('DEBUG:RegistrationFun:autoGenerateFOVmask: Would you like to see the final, auto-generated FOV mask?');
                if(showDebug)
                    msg = sprintf('DEBUG:RegistrationFun:autoGenerateFOVmask: Mask image generated from smoothed\nbackground, cutoff at %f, which is %f * the mean intensity of %f', illuminatedRegionCutoff, GlobalVars.maskIntensityCutoff, meanIntensity);
                    PlotFun.plotClusterImage(maskImage, msg);
                end
                drawnow;
            end
            %### END DEBUG ###
            
            %write the mask image to a file
            RegistrationFun.saveFOVmaskImage(maskImage, maskFilename);
        end
        
        
        function subset = identifyFilterSubsetsForRegistration(cDat, filterSubsets)
            if(isempty(filterSubsets))
                subset = true(cDat.numClusters,1); %use all clusters for registration if no filters were specified
            else
                subset = false(cDat.numClusters,1);
                for i = 1:size(filterSubsets,2)
                    if(cDat.isValidFilterName(filterSubsets{i}))
                        subset = subset | cDat.getFilteredSubset(filterSubsets{i});
                    else
                        errMsg = ['"' filterSubsets{i} '" is not a valid filter subset that is found in the sequence file'];
                        warning('RegistrationFun:identifyFilterSubsetsForRegistration',errMsg);
                    end
                end
            end
            if(size(find(subset),1) < GlobalVars.minNumClustersForRegistration) %if the filters specified did not contain enough clusters for registration
                error('RegistrationFun:identifyFilterSubsetsForRegistration','This sequence file with specified filtered subsets did not contain enough clusters for global registration');
            end
        end
        
        
        function saveFOVmaskImage(maskImage, maskImageFilename)
            imwrite(uint8(maskImage), maskImageFilename);
        end
        
        function maskImage = loadFOVmaskImage(maskImageFilename)
            maskImage = uint8(imread(maskImageFilename));
        end % identifyFilterSubsetsForRegistration
        
        % Register two images to get offset between them. Allows for global
        % registration by registering one reference image, and then finding
        % offset between reference image and actual image to analyse.
        function [offsetY, offsetX, registrationSuccess, correlationPeakHeight] = RegisterTileByImage(rawImage, refImage)
           

            % Normalized crosscorrelation of the subtile image versus the synthetic image
            corrMatrix = normxcorr2(rawImage, refImage);

            %fit the correlation matrix to a gaussian to find the peak
            %correlation to sub-pixel resolution.  The correlation peak height is
            %used as a measure of registration success
            maxOffset = [Inf Inf];
            [registrationSuccess yPeak xPeak correlationPeakHeight] = RegistrationFun.fitCorrelationGaussian(corrMatrix, maxOffset);

            if(registrationSuccess) %if the cross-correlation for this tile was successful
                corrOffset = [(yPeak - size(rawImage,1)) (xPeak - size(rawImage,2))]; %the offset of the correlation peak relative to the subtile image size

                % sum the relative offsets of the position of the central subimage (synthetic data) registered onto the actual image
                offsetY = corrOffset(1) - 1;
                offsetX = corrOffset(2) - 1;


            else
                %do not adjust cluster positions if registration failed
                offsetY = 0;
                offsetX = 0;

            end
        end % END RegisterTileByImage
               
        
        
     end % END static methods
end % END class RegistrationFun
