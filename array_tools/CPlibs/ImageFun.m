% Utility functions for manipulating and subdividing images, as well as
% creating synthetic images from data

classdef ImageFun
    methods (Static)

        %create a synthetic image from a list of cluster centers
        function synthImagePoints = makeSynthImagePoints(clusters_y, clusters_x, size_y, size_x,intensityValue)
            synthImagePoints = zeros(size_y, size_x); %initialize a matrix for the synthetic image and set the background to black (0)
            idx = sub2ind(size(synthImagePoints), round(clusters_y - min(clusters_y) + 1), round(clusters_x - min(clusters_x) + 1)); %determine the pixel indices on the synthetic image, when represented as a linear array, that correspond to cluster centers from the fit
            synthImagePoints(idx) = intensityValue;
        end
        
        %create a synthetic image from a list of 2D gaussians
        function synthImageGauss = makeSynthImageGauss(cDat,subIdx,imageBounds)
            n = size(subIdx,1); %number of gaussians to be plotted

            leftX = imageBounds(1);
            topY = imageBounds(2);
            xMaxGlobal = imageBounds(3)+1;
            yMaxGlobal = imageBounds(4)+1;
            
            synthImageGauss = zeros(yMaxGlobal, xMaxGlobal); %initialize the image array

            % get parameters for each gaussian
            amp = cDat.fitAmplitudes(subIdx);
            yCent = cDat.fitY(subIdx)-topY+1;
            xCent = cDat.fitX(subIdx)-leftX+1;
            sigma = cDat.fitSigmas(subIdx);

            %determine the local window in which each individual gaussian will be computed.
            %the tails of a normal distribution go to infinity, but for
            %computational efficiency, we compute each gaussian in a local
            %window, the size of which is chosen so that no signal is lost 
            %that is higher than the global "AmplitudeCutoffForlocalWindow" factor
            %below the peak amplitude

            %The localPixelCutoff is how far from the center this amplitude
            %value (on the y axis) is in pixels (on the x axis)
            if(n > 0)
                localPixelCutoff = calcLocalCutoffGaussian(GlobalVars.AmplitudeCutoffForLocalWindow, max(yMaxGlobal,xMaxGlobal), max(sigma));
                %create the local image matrix
                yMaxLocal = localPixelCutoff*2 + 1;
                xMaxLocal = localPixelCutoff*2 + 1;

                % loop through gaussians
                for i = 1:n %for each gaussian

                    %create a coordinate grid in x and y for distance calculations.  
                    %Because gaussian centers are in pixel space, the difference between
                    %the coordinate of a given pixel and the coordinate of the gaussian 
                    %center is the distance between them. 
                    [xCoordMatrix,yCoordMatrix] = meshgrid(1:xMaxLocal,1:yMaxLocal);
                    yCoordMatrix = yCoordMatrix - 1 + round(yCent(i)) - localPixelCutoff;
                    xCoordMatrix = xCoordMatrix - 1 + round(xCent(i)) - localPixelCutoff;

                    %create the image of the current gaussian on the local image matrix.
                    localImageGauss = amp(i)*exp(-((((yCoordMatrix-yCent(i)).^2) + ((xCoordMatrix-xCent(i)).^2))/(2*(sigma(i)^2)))); 

                    %add the local image to the corresponding subset region of the full image
                    globalYrange = (round(yCent(i))-localPixelCutoff):(round(yCent(i))+localPixelCutoff);
                    globalXrange = (round(xCent(i))-localPixelCutoff):(round(xCent(i))+localPixelCutoff);
                    yIndices = find((globalYrange > 0) & (globalYrange <= yMaxGlobal));
                    xIndices = find((globalXrange > 0) & (globalXrange <= xMaxGlobal));
                    globalYrange = globalYrange(yIndices);
                    globalXrange = globalXrange(xIndices);
                    localYrange = 1:yMaxLocal;
                    localXrange = 1:xMaxLocal;
                    localYrange = localYrange(yIndices);
                    localXrange = localXrange(xIndices);
                    synthImageGauss(globalYrange,globalXrange) = synthImageGauss(globalYrange,globalXrange) + localImageGauss(localYrange,localXrange);
                end
                %do not add background, it will be added externally to this function
            else
                %if no clusters are present in this (sub)image, return a black (0) image
                synthImageGauss = 0;
            end
        end
        
        %get the histogram of pixel values from an image
        function [binHeights binCenters] = getImageHistogram(image, numBins)
            if(~exist('numBins', 'var'))
                numBins = max(image(:)) - min(image(:));
                if(numBins < 10) %make sure we at least have some minimum number of bins
                    numBins = 10;
                end
            end
            
            [binHeights,binCenters] = hist(image(:),numBins); %get histogram

            %pixel intensity values jump between discreet values due to bit
            %conversion and skip bins, leaving uninformative "zero" bins 
            %interspersed between informative bins.  Remove all zero bins
            populatedIndices = find(binHeights > 0);
            binHeights = binHeights(populatedIndices);
            binCenters = binCenters(populatedIndices);
        end %END getImageHistogram
               
        %save cluster image as a 16bit uncompressed tiff
        function saveClusterImage(image, filename)
            imwrite(uint16(image),filename,'Compression','none','tif');
        end
        
        %get the background values from a cluster image, output pixel by
        %pixel into a new background image
        function backgroundImage = getBackgroundFromClusterImage(rawImage, maxClusterSize)
            if(~exist('maxClusterSize', 'var'))
                maxClusterSize = 10; % default size in cluster (this parameter has been optimized for images taken on the GAIIx or other imaging station with a coolsnap K4 camera and ~20x objective)
            end
            backgroundImage = imopen(rawImage,strel('disk',maxClusterSize)); %identify background by morphological opening of a disk shape larger than the maximum feature size
        
%             localDebug = false;
%             if(GlobalVars.debug)
%                 response = questdlg('DEBUG:ImageFun:getBackgroundFromClusterImage: Would you like to see debugging figures for image background subtraction?.','DEBUG','yes','no','no');
%                 waitfor(response);
%                 if(strcmp(response,'yes'))
%                     localDebug = true;
%                 end
%             end
%             
%             if(localDebug)
%                 PlotFun.plotClusterImage(rawImage, 'DEBUG:ImageFun:getBackgroundFromClusterImage: Raw Image');
%                 [binHeights binCenters] = ImageFun.getImageHistogram(rawImage);
%                 PlotFun.plotImageHistogram(binHeights,binCenters,'DEBUG:ImageFun:getBackgroundFromClusterImage: Raw Image');
%                 PlotFun.plotImageIntensitySurface(rawImage, 'DEBUG:ImageFun:getBackgroundFromClusterImage: Raw Image');
% 
%                 PlotFun.plotClusterImage(backgroundImage, 'DEBUG:ImageFun:getBackgroundFromClusterImage: Background Image');
%                 [binHeights binCenters] = ImageFun.getImageHistogram(backgroundImage);
%                 PlotFun.plotImageHistogram(binHeights,binCenters,'DEBUG:ImageFun:getBackgroundFromClusterImage: Background Image');
%                 PlotFun.plotImageIntensitySurface(backgroundImage, 'DEBUG:ImageFun:getBackgroundFromClusterImage: Background Image');
% 
%                 backgroundSubtractedImage = rawImage - backgroundImage;
%                 PlotFun.plotClusterImage(backgroundSubtractedImage, 'DEBUG:ImageFun:getBackgroundFromClusterImage: Background Subtracted Image');
%                 [binHeights binCenters] = ImageFun.getImageHistogram(backgroundSubtractedImage);
%                 PlotFun.plotImageHistogram(binHeights,binCenters,'DEBUG:ImageFun:getBackgroundFromClusterImage: Background Subtracted Image');
%                 PlotFun.plotImageIntensitySurface(backgroundSubtractedImage, 'DEBUG:ImageFun:getBackgroundFromClusterImage: Background Subtracted Image');
%             end

        end %END getBackgroundFromClusterImage
         
        %returns the coords of a central subtile of an image
        function [centralSubtileCoords offset] = getCentralSubtileCoords(imageSize, subtileHeight, subtileWidth)
            if(~exist('subtileHeight','var'))
                subtileHeight = GlobalVars.globalRegistrationSubtileHeight;
            end
            if(~exist('subtileWidth','var'))
                subtileWidth = GlobalVars.globalRegistrationSubtileWidth;
            end
            
            %the subtile cannot be larger than the image.  Check bounds and
            %trim to the full dimension of the image, if necessary
            if(subtileHeight > imageSize(1))
                subtileHeight = imageSize(1);
            end
            if(subtileHeight < 2)
                subtileHeight = 2;
            end
            if(subtileWidth > imageSize(2))
                subtileWidth = imageSize(2);
            end
            if(subtileWidth < 2)
                subtileWidth = 2;
            end
            
            %calculate the position of the upper left pixel
            minY = 1 + floor((imageSize(1)/2)-(subtileHeight/2));
            minX = 1 + floor((imageSize(2)/2)-(subtileWidth/2));
            
            %return [minX minY maxX maxY]
            centralSubtileCoords = [minX minY subtileWidth-1 subtileHeight-1];
            
            %return offset of central tile from original image
            offset = [minY-1 minX-1];
        end
        
        function [centralSubtile offset] = getCentralSubtileImage(image, subtileHeight, subtileWidth)
            [subtileCoords offset] = ImageFun.getCentralSubtileCoords(size(image), subtileHeight, subtileWidth);%get central subtile coords
            centralSubtile = imcrop(image,subtileCoords);%crop image to subtile
        end
        
        % Divide an image as closely as possible into a subtiles of specified dimensions
        function [innerSubtileMap outerSubtileMap, numSubtiles] = getOverlappingSubtileMapCoords(imageSize, innerWidthY, innerWidthX, overlapY, overlapX)
            imageSizeY = imageSize(1);
            imageSizeX = imageSize(2);
            
            %take dimensions of the are inside an outer border of the overlap width
            centralImageSizeY = imageSizeY - (2*overlapY);
            centralImageSizeX = imageSizeX - (2*overlapX);
            
            if((centralImageSizeY < 0) || (centralImageSizeX < 0))
                error('ImageFun:getOverlappingSubtileMapCoords','image too small to be tiled with this overlap');
            end
            if(innerWidthY > imageSizeY)
                innerWidthY = imageSizeY;
            end
            if(innerWidthX > imageSizeX)
                innerWidthX = imageSizeX;
            end

            numSubTilesY = floor(centralImageSizeY/innerWidthY);
            numSubTilesX = floor(centralImageSizeX/innerWidthX);
            
            innerSubtileMap = ImageFun.getSubtileMapCoords([centralImageSizeY centralImageSizeX], numSubTilesY, numSubTilesX);
            numSubtiles = size(innerSubtileMap,1);
            repmat([overlapX overlapY 0 0],numSubtiles,1);
            innerSubtileMap = innerSubtileMap + repmat([overlapX overlapY 0 0],size(innerSubtileMap,1),1); %adjust position for overlap border
            outerSubtileMap = innerSubtileMap + repmat([-overlapX -overlapY (2*overlapX) (2*overlapY)],size(innerSubtileMap,1),1); %expand the inner subtiles by the overlap to give the outer subtile map
        end
        
        % Divide an image as evenly as possible into a specified number of subtiles
        function subtileMap = getSubtileMapCoords(imageSize, numSubTilesY, numSubTilesX)
            imageSizeY = imageSize(1);
            imageSizeX = imageSize(2);

            if((imageSizeY < 0) || (imageSizeX < 0))
                error('ImageFun:getSubtileMapCoords','invalid image with height or width equal to zero');
            end
            if((imageSizeY < numSubTilesY) || (imageSizeX < numSubTilesX))
                error('ImageFun:getSubtileMapCoords','subtile map cannot contain more subdivisions in a dimension than the number of pixels');
            end

            %divide to get unadjusted widths
            widthY = floor(imageSizeY / numSubTilesY);
            widthX = floor(imageSizeX / numSubTilesX);

            %get remainders left over
            modY = mod(imageSizeY, numSubTilesY);
            modX = mod(imageSizeX, numSubTilesX);

            %divide the remainder across a subset of the tiles
            subtileBoundsY = 1:widthY:imageSizeY;
            if(modY ~= 0)
                pixelAdjustment = 1:modY-1;
                subtileBoundsY(end-modY+2:end) = subtileBoundsY(end-modY+2:end) + pixelAdjustment;
            else
                subtileBoundsY = [subtileBoundsY imageSizeY];
            end

            subtileBoundsX = 1:widthX:imageSizeX;
            if(modX ~= 0)
                pixelAdjustment = 1:modX-1;
                subtileBoundsX(end-modX+2:end) = subtileBoundsX(end-modX+2:end) + pixelAdjustment;
            else
                subtileBoundsX = [subtileBoundsX imageSizeX];
            end

            %create a array and populate with subtile map data (each entry a
            %vector of length 4: [leftPixel topPixel widthX heightY]
            subtileMap = zeros(numSubTilesY*numSubTilesX,4); %
            subtileIdx = 1;
            for i = 1:numSubTilesY
                topPixel = subtileBoundsY(i);
                heightY = subtileBoundsY(i+1) - subtileBoundsY(i);
                for j = 1:numSubTilesX
                    leftPixel = subtileBoundsX(j);
                    widthX = subtileBoundsX(j+1) - subtileBoundsX(j);
                    subtileMap(subtileIdx,:) = [leftPixel topPixel widthX heightY];
                    subtileIdx = subtileIdx + 1;
                end
            end
        end     
        
        
        
        %set a border around the edge of a matrix (e.g. an image) to a given pixel value
        function image = makeBorder(image, borderWidth, pixelValue)
            if(~exist('pixelValue','var') || isempty(pixelValue))
                pixelValue = 0; %default 0 (black)
            end
            image(1:borderWidth,:) = pixelValue;
            image(:,1:borderWidth) = pixelValue;
            image(end-borderWidth+1:end,:) = pixelValue;
            image(:,end-borderWidth+1:end) = pixelValue;
        end
        
    end % END STATIC METHODS
    
end % END ImageFun






%%%%%%%%%%% LOCAL METHODS %%%%%%%%%%%%%%

%given some fraction of the peak gaussian amplitude (in y), this function will
%give the distance in x away from the center
function localPixelCutoff = calcLocalCutoffGaussian(AmplitudeCutoffForlocalWindow, maxPixelValue, maxSigma)
    minPixelCutoff = 5; %minimum allowed local pixel window

    for i = maxPixelValue:-1:1
        currHeight = exp(-1*i^2/(2*(maxSigma^2)));
        if(currHeight ~= 0)
            if((1/currHeight) < AmplitudeCutoffForlocalWindow)
                break;
            elseif(i==1)
                break;
            end
        end
    end
    localPixelCutoff = i + 1;
    if(localPixelCutoff < minPixelCutoff)
        localPixelCutoff = minPixelCutoff;
    end
        
end