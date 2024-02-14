% Class for plotting images, data registrations, and fits to clustered
% flowcell images (e.g. Illumina)

classdef PlotFun
    
    properties (Constant)
        colorMap = containers.Map( {'cyan',...
                                    'red',...
                                    'blue',...
                                    'green',...
                                    'magenta',...
                                    'yellow',...
                                    'orange'},...
                                   {[0 1 1],...
                                    [1 0 0],...
                                    [0 0 1],...
                                    [0 1 0],...
                                    [1 0 1],...
                                    [1 1 0],...
                                    [1 0.647059 0]});
    end % END static properties
    
    methods (Static)
         
         %plot an actual (sub)image, and overlay a list of points e.g. to manually
         %display the quality of registration
         function plotRegistration(image, dataYX, bg, amp, figureTitle, figureHandle, axisHandle)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('bg', 'var') || isempty(bg))
                bg = min(image(:));
            end
            if(~exist('amp', 'var') || isempty(amp))
                amp = max(image(:))/GlobalVars.goldenMean;
                if(amp < bg)
                    amp = max(image(:));
                end
            end
            if(~exist('figureTitle', 'var'))
                figureTitle = [];
            end
            %plot the image
            [figureHandle axisHandle] = PlotFun.plotClusterImage(image, figureTitle, bg, amp, figureHandle, axisHandle);
            hold on
            %plot the registered data points
            plot(axisHandle,dataYX(:,2),dataYX(:,1),'+','linestyle','none','Color','red','MarkerSize',2);
         end % END plotRegistration
         
         
         function colorName = getNextColor(index)
            colorList{2}  = PlotFun.colorMap('cyan');
            colorList{3}  = PlotFun.colorMap('red');
            colorList{4}  = PlotFun.colorMap('blue');
            colorList{5}  = PlotFun.colorMap('green');
            colorList{6}  = PlotFun.colorMap('magenta');
            colorList{7}  = PlotFun.colorMap('yellow');
            colorList{1}  = PlotFun.colorMap('orange');
            colorName = colorList{mod(index,7)+1};
         end
         
         function plotGlobalRegistrationByFilter(image, cDat, bg, amp, figureTitle, filterColormap, figureHandle, axisHandle)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('bg', 'var') || isempty(bg))
                bg = min(image(:));
            end
            if(~exist('amp', 'var') || isempty(amp))
                amp = max(image(:))/GlobalVars.goldenMean;
                if(amp < bg)
                    amp = max(image(:));
                end
            end
            if(~exist('figureTitle', 'var'))
                figureTitle = [];
            end
            if(~exist('filterColormap','var'))
                filterColormap = containers.Map();
            end
            
            %plot the image
            [figureHandle axisHandle] = PlotFun.plotClusterImage(image, figureTitle, bg, amp, figureHandle, axisHandle);
            hold on

            %plot registered data points for each filter
            filterIDs = keys(cDat.filterSubsets); %get array of filter IDs
            legendIndex = 1;
            for i = 1:length(filterIDs)
                currY = cDat.adjY(~cDat.excised & cDat.filterSubsets(filterIDs{i}));
                currX = cDat.adjX(~cDat.excised & cDat.filterSubsets(filterIDs{i}));
                if isempty(keys(filterColormap)) % if a filterColormap wasn't defined, assign colors using getNextColor function
                    currColor = PlotFun.getNextColor(i);
                else % if a filterColormap is defined, assign the color based on what is specified in filterColormap
                    if isKey(filterColormap, filterIDs{i})
                        currColor = PlotFun.colorMap(filterColormap(filterIDs{i}));
                    else % this filter does not have a specified color in the colormap, so look for the special ELSE clause
                        if isKey(filterColormap, 'ELSE')
                            currColor = PlotFun.colorMap(filterColormap('ELSE'));
                        else
                            currColor = PlotFun.colorMap('blue'); % if 'ELSE' isn't defined in the filterColormap, just use the color 'blue'
                        end
                    end
                end
                if(~isempty(currY))
                    plot(axisHandle,currX,currY,'+','Linestyle','None','MarkerSize',3,'Color',currColor,'MarkerFaceColor',currColor);
                    legendTitles{legendIndex} = strtrim(strrep(filterIDs{i},'_',' ')); %change underscores (which are used in subscript syntax) to spaces
                    legendIndex = legendIndex + 1;
                end
            end
            if(legendIndex > 1) % show the legend
                legend(axisHandle,legendTitles);
            end
         end
         
         
         %plot a synthetic image from a 2D fits
         function plotFitImage(fitImage, dataYX, bg, amp, figureTitle, figureHandle, axisHandle)
             if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('bg', 'var') || isempty(bg))
                bg = min(fitImage(:));
            end
            if(~exist('amp', 'var') || isempty(amp))
                amp = max(fitImage(:))/GlobalVars.goldenMean;
                if(amp < bg)
                    amp = max(fitImage(:));
                end
            end
            
            %plot the synthetic fit image
            [figureHandle axisHandle] = PlotFun.plotClusterImage(fitImage, figureTitle, bg, amp, figureHandle, axisHandle);
            hold on
            
            %plot the qseq data points overlaid on the fit image
            plot(axisHandle,dataYX(:,2),dataYX(:,1),'+','linestyle','none','Color','red','MarkerSize',2);
         end % END plotFitImage
         
         function figureHandle = saveThreePanelFitWithoutGui(fullImage ,cDat, bgImage, subimageBounds, dataBounds, qcFilename, globalOffset_y, globalOffset_x)
            figureHandle = figure('visible','off');
            
            if(~exist('globalOffset_x','var') || isempty(globalOffset_x || ~exist('globalOffset_y','var') || isempty(globalOffset_y)))
                figureHandle = PlotFun.plotThreePanelFit(fullImage, cDat, bgImage, subimageBounds, dataBounds, qcFilename, [], figureHandle);
            else
                figureHandle = PlotFun.plotThreePanelFit(fullImage, cDat, bgImage, subimageBounds, dataBounds, qcFilename, [], figureHandle, globalOffset_y, globalOffset_x);
            end
            
            close(figureHandle);
         end
                 
         function figureHandle = plotThreePanelFit(fullImage, cDat, bgImage, subimageBounds, dataBounds, qcFilename, residualBorder, figureHandle, globalOffset_y, globalOffset_x)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure();
            end
            if(~exist('bgImage','var') || isempty(bgImage))
                bgImage = zeros(size(fullImage)); %background to be added into the synthetic image
            else
                if(isscalar(bgImage))
                    %if the background was passed in as a scalar, make a
                    %bgImage matrix of the appropriate size
                    bgImage = bgImage*ones(size(fullImage)); 
                elseif(~isequal(size(bgImage),size(fullImage)))
                    error('PlotFun:plotThreePanelFit','background image must be the same size as the full image.');
                end
            end          
            if(~exist('subimageBounds','var') || isempty(subimageBounds))
                subimageBounds = ImageFun.getCentralSubtileCoords(size(fullImage), 300, 300);
            end
            
            imageLeftX = subimageBounds(1);
            imageTopY = subimageBounds(2);
            
            if(~exist('dataBounds','var') || isempty(dataBounds))
                dataBounds = subimageBounds;
            end
            if(~exist('residualBorder','var') || isempty(residualBorder))
                residualBorder = 0;
            end
            
            if(~exist('globalOffset_x','var') || isempty(globalOffset_x || ~exist('globalOffset_y','var') || isempty(globalOffset_y)))
                printOffsets = false;
            else
                printOffsets = true;
            end
            
            dataLeftX = dataBounds(1);
            dataRightX = dataBounds(1) + dataBounds(3);
            dataTopY = dataBounds(2);
            dataBottomY = dataBounds(2) + dataBounds(4);
            
            %set the figure size
            set(figureHandle, 'Position', [100, 100, 1600, 800]);
            
            %crop image to the subimage region
            subImage = imcrop(fullImage,subimageBounds);
            bgImage = imcrop(bgImage,subimageBounds);
           
            %only use fit data from (unexcised) points within the subimage region
            subIdx = find(~cDat.excised & cDat.adjY>=dataTopY & cDat.adjY<dataBottomY & cDat.adjX>=dataLeftX & cDat.adjX<dataRightX); %get indices that fall within the subimage region

            %get the image scaling
            subImageBg = min(subImage(:));
            subImageAmp = max(subImage(:))/GlobalVars.goldenMean;
            if(subImageAmp < subImageBg)
                subImageAmp = max(subImage(:));
            end
            
            % plot the fit cluster centers on the acutal subimage in a subplot window (left)
            axisHandleLeft = subplot(1,3,1, 'Parent', figureHandle); %left
            dataYX = [cDat.adjY(subIdx)-imageTopY+1 cDat.adjX(subIdx)-imageLeftX+1];
            PlotFun.plotRegistration(subImage, dataYX, subImageBg, subImageAmp, 'actual image (registered)', figureHandle, axisHandleLeft);
            colorbarHandleLeft = colorbar('peer',axisHandleLeft,'location','SouthOutside');
            cbfreeze(colorbarHandleLeft);
            freezeColors(axisHandleLeft);
            hold on
            if(printOffsets)
                text(2,-50,{strcat('xOffSet: ',num2str(globalOffset_x)),strcat('yOffSet: ',num2str(globalOffset_y))}, 'color','red','FontSize',12,'FontName','Arial','FontWeight','bold','VerticalAlignment','Top','Interpreter','none');
            end
   
            % create the synthetic fit image reconstructed from 2D gaussians in a subplot window (middle)
            fitImage = bgImage + ImageFun.makeSynthImageGauss(cDat,subIdx,subimageBounds); %generate the synthetic image
            axisHandleMiddle = subplot(1,3,2, 'Parent', figureHandle); %middle
            PlotFun.plotFitImage(fitImage, dataYX, subImageBg, subImageAmp, 'fit image', figureHandle, axisHandleMiddle);
            colorbarHandleMiddle = colorbar('peer',axisHandleMiddle,'location','SouthOutside');
            cbfreeze(colorbarHandleMiddle);
            freezeColors(axisHandleMiddle);
            
            % plot the residuals from the fit (image - fit) in a subplot window (right)
            axisHandleRight = subplot(1,3,3, 'Parent', figureHandle); %right
            residualImage = double(subImage)-round(fitImage);
            residualImage = ImageFun.makeBorder(residualImage, residualBorder);
            imagesc(residualImage,'Parent', axisHandleRight,[-max(fitImage(:))/GlobalVars.goldenMean max(fitImage(:))/GlobalVars.goldenMean]);
            axis(axisHandleRight, 'square');
            colormap(axisHandleRight, 'jet');
            colorbarHandleRight = colorbar('peer',axisHandleRight,'location','SouthOutside');
            cbfreeze(colorbarHandleRight);
            freezeColors(axisHandleRight);
            set(axisHandleRight,'YTick',[]);
            set(axisHandleRight,'XTick',[]);
            title(axisHandleRight,'residual (actual - fit)');

            % if a filename is provided, print the figure to a file
            if(exist('qcFilename','var') && ~isempty(qcFilename))
                if(ishandle(figureHandle))
                    set(figureHandle,'PaperPositionMode','auto');                   
                    print(figureHandle,'-depsc','-painters','-r300', qcFilename); %export the figure as a .eps
                    disp(['Three-panel fit QC figure saved to: ' qcFilename]);
                end
            end
         end

         
        %plot image intensity values on a 3D surface
        function [figureHandle axisHandle surfaceplotHandle] = plotImageIntensitySurface(image, axisTitle, interval, figureHandle, axisHandle)
            %plotting an entire image is often too many points
            %"interval" allows only every Nth (e.g. 10th) pixel value to be
            %plotted, under the assumption that local regions have similar intensities
            if(~exist('interval', 'var'))
                %create a 2D x,y grid to map cluster shift data (z axis).
                %The grid will have a discreet point every "interval'th" pixel
                interval = 1 + round(sqrt(numel(image)/1e5));
            end
            if(~exist('axisTitle', 'var'))
                axisTitle = '';
            end
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            surfaceplotHandle = surf(axisHandle, double(image(1:interval:end,1:interval:end)));
            zlim(axisHandle,[min(image(:)) max(image(:))]);
            imageSize = size(image);
            ylim(axisHandle,[0 ceil(imageSize(1)/interval)]);
            xlim(axisHandle,[0 ceil(imageSize(2)/interval)]);
            set(axisHandle,'ydir','reverse');
            %scale the axes to account for the sampling interval
            oldYticks = get(axisHandle,'ytick');
            set(axisHandle,'yticklabel',num2str(round(oldYticks'*interval)));
            oldXticks = get(axisHandle,'xtick');
            set(axisHandle,'xticklabel',num2str(round(oldXticks'*interval)));

            title(axisHandle, axisTitle);
        end %END plotImageIntensitySurface
        
        %plot grayscale 16bit cluster image in a figure window
        function [figureHandle axisHandle imageHandle] = plotClusterImage(image, axisTitle, minVal, maxVal, figureHandle, axisHandle)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('axisTitle', 'var'))
                axisTitle = '';
            end
            if(~exist('minVal', 'var') || isempty(minVal))
                minVal = min(image(:)); 
            end
            if(~exist('maxVal', 'var') || isempty(maxVal))
                maxVal = max(image(:))/GlobalVars.goldenMean;
                if(maxVal < minVal)
                    maxVal = max(image(:));
                end
            end
            imageHandle = imagesc(image,'Parent',axisHandle,[minVal maxVal]);
            set(axisHandle,'YTick',[]);
            set(axisHandle,'XTick',[]);
            axis(axisHandle,'square');
            colormap(axisHandle,'gray');
            title(axisHandle, axisTitle);
        end %END plotImage
         
        %plot an image histogram
        function [figureHandle axisHandle areaseriesHandle] = plotImageHistogram(binHeights, binCenters, axisTitle, minAxis, maxAxis, figureHandle, axisHandle)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('axisTitle', 'var'))
                axisTitle = '';
            end
            if(~exist('minAxis', 'var'))
                minAxis = min(binCenters);
            end
            if(~exist('maxAxis', 'var'))
                maxAxis = max(binCenters);
            end
            areaseriesHandle = area(axisHandle,binCenters,binHeights);
            axis(axisHandle,[minAxis maxAxis 0 max(binHeights)]);
            colormap(axisHandle,'gray');
            ylabel(axisHandle,'# pixels (histogram)');
            xlabel(axisHandle,'intensity value');
            title(axisHandle, axisTitle);
        end

        function plotFitCorrelationMatrix(corrMatrix,minY,minX,params,figureTitle,figureHandle,axisHandle,interval)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('interval', 'var'))
                interval = 1 + round(sqrt(numel(corrMatrix)/1e5));
            end
            if(~exist('figureTitle', 'var'))
                figureTitle = '';
            end
            if(~exist('minY', 'var'))
                minY = 1;
            end
            if(~exist('minX', 'var'))
                minX = 1;
            end
            
            %plot correlation matrix (in color)
            [figureHandle axisHandle] = PlotFun.plotCorrelationMatrix(corrMatrix,minY,minX,figureTitle,figureHandle,axisHandle,interval);
            freezeColors(axisHandle);
            hold on;
            
            %plot fit surface (in gray)
            sizeCorrMatrix = size(corrMatrix);
            [xq,yq] = meshgrid(1:interval:sizeCorrMatrix(2),1:interval:sizeCorrMatrix(1)); %make 2d mesh grid arrays (where every matrix position contains its index value)
            yxq(:,1) = yq(:) + minY - 1;
            yxq(:,2) = xq(:) + minX - 1;

            zq = FitFun.gaussian_2D(params,yxq);
            zqMat = vec2mat(zq,sizeCorrMatrix(2));
            colormap(axisHandle,'gray');
            surfaceplotHandle = mesh(axisHandle,yq,xq,zqMat); %plot the fit on a 3D mesh plot
            freezeColors(axisHandle);
        end
        
        function [figureHandle axisHandle surfaceplotHandle] = plotCorrelationMatrix(corrMatrix,minY,minX,figureTitle,figureHandle,axisHandle,interval)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('interval', 'var'))
                interval = 1 + round(sqrt(numel(corrMatrix)/1e5));
            end
            if(~exist('figureTitle', 'var'))
                figureTitle = '';
            end
            if(~exist('minY', 'var'))
                minY = 1;
            end
            if(~exist('minX', 'var'))
                minX = 1;
            end
            
            surfaceplotHandle = surf(axisHandle, corrMatrix(1:interval:end,1:interval:end));
            zlim(axisHandle,[min(corrMatrix(:)) max(corrMatrix(:))]);
            corrMatrixSize = size(corrMatrix);
            ylim(axisHandle,[0 ceil(corrMatrixSize(1)/interval)]);
            xlim(axisHandle,[0 ceil(corrMatrixSize(2)/interval)]);        
            set(axisHandle,'ydir','reverse');
            
            %scale the axes to account for the sampling interval and offset
            %if this is a subregion of a larger correlation matrix
            oldYticks = get(axisHandle,'ytick');
            set(axisHandle,'yticklabel',num2str(minY - 1 + round(oldYticks'*interval)));
            set(axisHandle, 'YTickMode','manual');
            oldXticks = get(axisHandle,'xtick');
            set(axisHandle,'xticklabel',num2str(minX - 1 + round(oldXticks'*interval)));
            set(axisHandle, 'XTickMode','manual');
            ylabel(axisHandle,'y');
            xlabel(axisHandle,'x');
            zlabel(axisHandle,'cross correlation');
            title(axisHandle, figureTitle);
            colormap(axisHandle,'jet');
            freezeColors(axisHandle)
            hold on;
        end
        
        function [figureHandle axisHandle surfaceplotHandle] = plotRegistrationOffsetFit(clusters_y, clusters_x, localOffsets, figureTitle, fitParams, figureHandle, interval)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('interval', 'var'))
                %create a 2D y,x grid to map cluster shift data (z axis).
                %The grid will have a discrete point every "interval'th" pixel
                interval = 1 + round(length(clusters_y)/2e4);
            end
            if(~exist('figureTitle', 'var'))
                figureTitle = '';
            end
            if(length(fitParams) > 2)
                params = fitParams(1:end-2);
                y0 = fitParams(end-1);
                x0 = fitParams(end);
            else
                error('PlotFun:adjust', 'not enough parameters passed in to define the registration offset maps');
            end
                     
            %plot offset data (in color)
            [figureHandle axisHandle] = PlotFun.plotRegistrationOffsets(clusters_y,clusters_x,localOffsets,figureTitle,figureHandle);

            %plot fit surface (in gray)
            [xq,yq] = meshgrid(min(clusters_x):interval:max(clusters_x),min(clusters_y):interval:max(clusters_y)); %make 2d mesh grid arrays (where every matrix position contains its index value)
            xyq = [yq(:) xq(:)];
            sizeXq = size(xq);
            zq = FitFun.quadraticSurface(params,xyq,y0,x0);
            zqMat = vec2mat(zq,sizeXq(1));
            colormap(axisHandle,'gray');
            surfaceplotHandle = mesh(axisHandle,xq,yq,zqMat'); %plot the fit on a 3D mesh plot
        end

        function [figureHandle axisHandle surfaceplotHandle] = plotRegistrationOffsets(clusters_y,clusters_x,localOffsets,figureTitle,figureHandle,axisHandle,interval)
            if(~exist('figureHandle', 'var') || ~ishandle(figureHandle))
                figureHandle = figure;
            end
            if(~exist('axisHandle', 'var') || ~ishandle(axisHandle))
                axisHandle = axes('Parent',figureHandle);
            end
            if(~exist('interval', 'var'))
                %create a 2D x,y grid to map cluster shift data (z axis).
                %The grid will have a discreet point every "interval'th" pixel
                interval = 1 + round(length(clusters_x)/2e4);
            end
            [xq,yq] = meshgrid(min(clusters_x):interval:max(clusters_x),min(clusters_y):interval:max(clusters_y)); %make 2d mesh grid arrays (where every matrix position contains its index value)
            vq = griddata(clusters_x,clusters_y,localOffsets,xq,yq); %interpolate the data onto the grid ("snap" the z values to the neares x,y coordinate position)

            surfaceplotHandle = mesh(axisHandle,xq,yq,vq); %plot the grid on a 3D mesh plot
            ylabel(axisHandle,'y');
            xlabel(axisHandle,'x');
            zlabel(axisHandle,'registration offset');
            title(axisHandle, figureTitle);
            colormap(axisHandle,'jet');
            freezeColors(axisHandle)
            hold on;
        end
        
        function [adjY adjX] = adjustRegistrationOffset(clusters_y,clusters_x,fitParamsY,fitParamsX)
            if((length(fitParamsY) > 2) && (length(fitParamsX) > 2))
                paramsY = fitParamsY(1:end-2);
                y0_y = fitParamsY(end-1);
                x0_y = fitParamsY(end);
                paramsX = fitParamsX(1:end-2);
                y0_x = fitParamsX(end-1);
                x0_x = fitParamsX(end);
            else
                error('PlotFun:adjustRegistrationOffset', 'not enough parameters passed in to define the registration offset maps');
            end
            
            offsetMapY = FitFun.quadraticSurface(paramsY,[clusters_y clusters_x],y0_y,x0_y);
            offsetMapX = FitFun.quadraticSurface(paramsX,[clusters_y clusters_x],y0_x,x0_x);
            
            if((y0_y == y0_x) && (x0_y == x0_x))
                y0 = y0_x;
                x0 = x0_x;
            else
                error('PlotFun:adjustRegistrationOffset','x and y registration offsets should have been calculated on the same data')
            end
            
            adjY = clusters_y + offsetMapY - y0;
            adjX = clusters_x + offsetMapX - x0;
        end
        
        function fitParams = fitRegistrationOffset(clusters_y, clusters_x, localOffsets)
            y0 = mean(clusters_y);
            x0 = mean(clusters_x);

            lb           = [-1e-3, -1e-3, -0.1, -0.1, -1e-3, min(localOffsets)];
            initialGuess = [    0,     0,    0,    0,    0, mean(localOffsets)];
            ub           = [ 1e-3,  1e-3,  0.1,  0.1,  1e-3, max(localOffsets)];
 
            options = optimset( 'TolX', 1e-12,... % TolX = convergence criteron on the change in parameters one iteration to the next
                            'TolFun', 1e-12,... % TolX = convergence criteron on the change in the value of the objective function one iteration to the next
                            'MaxFunEvals', 10000,... %max number of times objective function can be evaluated
                            'MaxIter', 1000,... %max number of iterations of the solver
                            'Display', 'off',... %display no output from optimization
                            'jacobian', 'on',... %use the jacobian
                            'DerivativeCheck', 'off',... %do not compare user-supplied jacobian (gradient of objective) to calculated finite-differencing derivatives
                            'LargeScale', 'on');

            fun = @(params,yx)FitFun.quadraticSurface(params,[clusters_y clusters_x],y0,x0);
            
            fitParams = lsqcurvefit(fun,initialGuess,[clusters_y clusters_x],localOffsets,lb,ub,options);
            
            fitParams = [fitParams y0 x0]; %add the y0 and x0 center points to the fit params
        end
 
     end %END STATIC METHODS
end


% ============================================================================ %
% Local functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% getCDataHandles -- get handles of all descendents with indexed CData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hout = getCDataHandles(h)
    % getCDataHandles  Find all objects with indexed CData

    %recursively descend object tree, finding objects with indexed CData
    % An exception: don't include children of objects that themselves have CData:
    %   for example, scattergroups are non-standard hggroups, with CData. Changing
    %   such a group's CData automatically changes the CData of its children,
    %   (as well as the children's handles), so there's no need to act on them.

    error(nargchk(1,1,nargin,'struct'))

    hout = [];
    if isempty(h),return;end

    ch = get(h,'children');
    for hh = ch'
        g = get(hh);
        if isfield(g,'CData'),     %does object have CData?
            %is it indexed/scaled?
            if ~isempty(g.CData) && isnumeric(g.CData) && size(g.CData,3)==1,
                hout = [hout; hh]; %#ok<AGROW> %yes, add to list
            end
        else %no CData, see if object has any interesting children
            hout = [hout; getCDataHandles(hh)]; %#ok<AGROW>
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% getParentAxes -- return handle of axes object to which a given object belongs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hAx = getParentAxes(h)
    % getParentAxes  Return enclosing axes of a given object (could be self)

    error(nargchk(1,1,nargin,'struct'))
    %object itself may be an axis
    if strcmp(get(h,'type'),'axes'),
        hAx = h;
        return
    end

    parent = get(h,'parent');
    if (strcmp(get(parent,'type'), 'axes')),
        hAx = parent;
    else
        hAx = getParentAxes(parent);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% checkArgs -- Validate input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, nancolor] = checkArgs(args)
    % checkArgs  Validate input arguments to freezeColors

    nargs = length(args);
    error(nargchk(0,3,nargs,'struct'))

    %grab handle from first argument if we have an odd number of arguments
    if mod(nargs,2),
        h = args{1};
        if ~ishandle(h),
            error('JRI:freezeColors:checkArgs:invalidHandle',...
                'The first argument must be a valid graphics handle (to an axis)')
        end
        % 4/2010 check if object to be frozen is a colorbar
        if strcmp(get(h,'Tag'),'Colorbar'),
            if ~exist('cbfreeze.m'),
                warning('JRI:freezeColors:checkArgs:cannotFreezeColorbar',...
                    ['You seem to be attempting to freeze a colorbar. This no longer'...
                    'works. Please read the help for freezeColors for the solution.'])
            else
                cbfreeze(h);
                return
            end
        end
        args{1} = [];
        nargs = nargs-1;
    else
        h = gca;
    end

    %set nancolor if that option was specified
    nancolor = [nan nan nan];
    if nargs == 2,
        if strcmpi(args{end-1},'nancolor'),
            nancolor = args{end};
            if ~all(size(nancolor)==[1 3]),
                error('JRI:freezeColors:checkArgs:badColorArgument',...
                    'nancolor must be [r g b] vector');
            end
            nancolor(nancolor>1) = 1; nancolor(nancolor<0) = 0;
        else
            error('JRI:freezeColors:checkArgs:unrecognizedOption',...
                'Unrecognized option (%s). Only ''nancolor'' is valid.',args{end-1})
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% freezeColors - lock the colors of a plot so multiple colormaps are
%%% possible between different subplots of the same figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freezeColors(varargin)
    % freezeColors  Lock colors of plot, enabling multiple colormaps per figure. (v2.3)
    %
    % Free for all uses, but please retain the following:
    %   Original Author:
    %   John Iversen, 2005-10
    %   john_iversen@post.harvard.edu

    appdatacode = 'JRI__freezeColorsData';

    [h, nancolor] = checkArgs(varargin);

    %gather all children with scaled or indexed CData
    cdatah = getCDataHandles(h);

    %current colormap
    cmap = colormap;
    nColors = size(cmap,1);
    cax = caxis;

    % convert object color indexes into colormap to true-color data using
    %  current colormap
    for hh = cdatah',
        g = get(hh);

        %preserve parent axis clim
        parentAx = getParentAxes(hh);
        originalClim = get(parentAx, 'clim');

        %   Note: Special handling of patches: For some reason, setting
        %   cdata on patches created by bar() yields an error,
        %   so instead we'll set facevertexcdata instead for patches.
        if ~strcmp(g.Type,'patch'),
            cdata = g.CData;
        else
            cdata = g.FaceVertexCData;
        end

        %get cdata mapping (most objects (except scattergroup) have it)
        if isfield(g,'CDataMapping'),
            scalemode = g.CDataMapping;
        else
            scalemode = 'scaled';
        end

        %save original indexed data for use with unfreezeColors
        siz = size(cdata);
        setappdata(hh, appdatacode, {cdata scalemode});

        %convert cdata to indexes into colormap
        if strcmp(scalemode,'scaled'),
            %4/19/06 JRI, Accommodate scaled display of integer cdata:
            %       in MATLAB, uint * double = uint, so must coerce cdata to double
            %       Thanks to O Yamashita for pointing this need out
            idx = ceil( (double(cdata) - cax(1)) / (cax(2)-cax(1)) * nColors);
        else %direct mapping
            idx = cdata;
            %10/8/09 in case direct data is non-int (e.g. image;freezeColors)
            % (Floor mimics how matlab converts data into colormap index.)
            % Thanks to D Armyr for the catch
            idx = floor(idx);
        end

        %clamp to [1, nColors]
        idx(idx<1) = 1;
        idx(idx>nColors) = nColors;

        %handle nans in idx
        nanmask = isnan(idx);
        idx(nanmask)=1; %temporarily replace w/ a valid colormap index

        %make true-color data--using current colormap
        realcolor = zeros(siz);
        for i = 1:3,
            c = cmap(idx,i);
            c = reshape(c,siz);
            c(nanmask) = nancolor(i); %restore Nan (or nancolor if specified)
            realcolor(:,:,i) = c;
        end

        %apply new true-color color data

        %true-color is not supported in painters renderer, so switch out of that
        if strcmp(get(gcf,'renderer'), 'painters'),
            set(gcf,'renderer','zbuffer');
        end

        %replace original CData with true-color data
        if ~strcmp(g.Type,'patch'),
            set(hh,'CData',realcolor);
        else
            set(hh,'faceVertexCData',permute(realcolor,[1 3 2]))
        end

        %restore clim (so colorbar will show correct limits)
        if ~isempty(parentAx),
            set(parentAx,'clim',originalClim)
        end

    end %loop on indexed-color objects
end % END freezeColors



function CBH = cbfreeze(varargin)
%CBFREEZE   Freezes the colormap of a colorbar.
%
%   SYNTAX:
%           cbfreeze
%           cbfreeze('off')
%           cbfreeze(H,...)
%     CBH = cbfreeze(...);
%
%   INPUT:
%     H     - Handles of colorbars to be freezed, or from figures to search
%             for them or from peer axes (see COLORBAR).
%             DEFAULT: gcf (freezes all colorbars from the current figure)
%     'off' - Unfreezes the colorbars, other options are:
%               'on'    Freezes
%               'un'    same as 'off'
%               'del'   Deletes the colormap(s).
%             DEFAULT: 'on' (of course)
%
%   OUTPUT (all optional):
%     CBH - Color bar handle(s).
%
%   DESCRIPTION:
%     MATLAB works with a unique COLORMAP by figure which is a big
%     limitation. Function FREEZECOLORS by John Iversen allows to use
%     different COLORMAPs in a single figure, but it fails freezing the
%     COLORBAR. This program handles this problem.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If no colorbar is found, one is created.
%     * The new frozen colorbar is an axes object and does not behaves
%       as normally colorbars when resizing the peer axes. Although, some
%       time the normal behavior is not that good.
%     * Besides, it does not have the 'Location' property anymore.
%     * But, it does acts normally: no ZOOM, no PAN, no ROTATE3D and no
%       mouse selectable.
%     * No need to say that CAXIS and COLORMAP must be defined before using
%       this function. Besides, the colorbar location. Anyway, 'off' or
%       'del' may help.
%     * The 'del' functionality may be used whether or not the colorbar(s)
%       is(are) froozen. The peer axes are resized back. Try: 
%        >> colorbar, cbfreeze del
%
%   EXAMPLE:
%     surf(peaks(30))
%     colormap jet
%     cbfreeze
%     colormap gray
%     title('What...?')
%
%   SEE ALSO:
%     COLORMAP, COLORBAR, CAXIS
%     and
%     FREEZECOLORS by John Iversen
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbfreeze.m
%   VERSION: 1.1 (Sep 02, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.1      Fixed BUG with image handle on MATLAB R2009a. Thanks to Sergio
%            Muniz. (Sep 02, 2009)

%   DISCLAIMER:
%   cbfreeze.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
cbappname = 'Frozen';         % Colorbar application data with fields:
                              % 'Location' from colorbar
                              % 'Position' from peer axes befor colorbar
                              % 'pax'      handle from peer axes.
axappname = 'FrozenColorbar'; % Peer axes application data with frozen
                              % colorbar handle.
 
% Set defaults:
S = 'on';                   Sopt = {'on','un','off','del'};
H = get(0,'CurrentFig');

% Check inputs:
if nargin==2 && (~isempty(varargin{1}) && all(ishandle(varargin{1})) && ...
  isempty(varargin{2}))
 
 % Check for CallBacks functionalities:
 % ------------------------------------
 
 varargin{1} = double(varargin{1});
 
 if strcmp(get(varargin{1},'BeingDelete'),'on') 
  % Working as DeletFcn:

  if (ishandle(get(varargin{1},'Parent')) && ...
      ~strcmpi(get(get(varargin{1},'Parent'),'BeingDeleted'),'on'))
    % The handle input is being deleted so do the colorbar:
    S = 'del'; 
    
   if ~isempty(getappdata(varargin{1},cbappname))
    % The frozen colorbar is being deleted:
    H = varargin{1};
   else
    % The peer axes is being deleted:
    H = ancestor(varargin{1},{'figure','uipanel'}); 
   end
   
  else
   % The figure is getting close:
   return
  end
  
 elseif (gca==varargin{1} && ...
                     gcf==ancestor(varargin{1},{'figure','uipanel'}))
  % Working as ButtonDownFcn:
  
  cbfreezedata = getappdata(varargin{1},cbappname);
  if ~isempty(cbfreezedata) 
   if ishandle(cbfreezedata.ax)
    % Turns the peer axes as current (ignores mouse click-over):
    set(gcf,'CurrentAxes',cbfreezedata.ax);
    return
   end
  else
   % Clears application data:
   rmappdata(varargin{1},cbappname) 
  end
  H = varargin{1};
 end
 
else
 
 % Checks for normal calling:
 % --------------------------
 
 % Looks for H:
 if nargin && ~isempty(varargin{1}) && all(ishandle(varargin{1}))
  H = varargin{1};
  varargin(1) = [];
 end

 % Looks for S:
 if ~isempty(varargin) && (isempty(varargin{1}) || ischar(varargin{1}))
  S = varargin{1};
 end
end

% Checks S:
if isempty(S)
 S = 'on';
end
S = lower(S);
iS = strmatch(S,Sopt);
if isempty(iS)
 error('CVARGAS:cbfreeze:IncorrectStringOption',...
  ['Unrecognized ''' S ''' argument.' ])
else
 S = Sopt{iS};
end

% Looks for CBH:
CBH = cbhandle(H);

if ~strcmp(S,'del') && isempty(CBH)
 % Creates a colorbar and peer axes:
 pax = gca;
 CBH = colorbar('peer',pax);
else
 pax = [];
end


% -------------------------------------------------------------------------
% MAIN 
% -------------------------------------------------------------------------
% Note: only CBH and S are necesary, but I use pax to avoid the use of the
%       "hidden" 'Axes' COLORBAR's property. Why... �?

% Saves current position:
fig = get(  0,'CurrentFigure');
cax = get(fig,'CurrentAxes');

% Works on every colorbar:
for icb = 1:length(CBH)
 
 % Colorbar axes handle:
 h  = double(CBH(icb));
 
 % This application data:
 cbfreezedata = getappdata(h,cbappname);
 
 % Gets peer axes:
 if ~isempty(cbfreezedata)
  pax = cbfreezedata.pax;
  if ~ishandle(pax) % just in case
   rmappdata(h,cbappname)
   continue
  end
 elseif isempty(pax) % not generated
  try
   pax = double(get(h,'Axes'));  % NEW feature in COLORBARs
  catch
   continue
  end
 end
 
 % Choose functionality:
 switch S
  
  case 'del'
   % Deletes:
   if ~isempty(cbfreezedata)
    % Returns axes to previous size:
    oldunits = get(pax,'Units');
    set(pax,'Units','Normalized');
    set(pax,'Position',cbfreezedata.Position)
    set(pax,'Units',oldunits)
    set(pax,'DeleteFcn','')
    if isappdata(pax,axappname)
     rmappdata(pax,axappname)
    end
   end
   if strcmp(get(h,'BeingDelete'),'off') 
    delete(h)
   end
   
  case {'un','off'}
   % Unfrozes:
   if ~isempty(cbfreezedata)
    delete(h);
    set(pax,'DeleteFcn','')
    if isappdata(pax,axappname)
     rmappdata(pax,axappname)
    end
    oldunits = get(pax,'Units');
    set(pax,'Units','Normalized')
    set(pax,'Position',cbfreezedata.Position)
    set(pax,'Units',oldunits)
    CBH(icb) = colorbar(...
     'peer'    ,pax,...
     'Location',cbfreezedata.Location);
   end
 
  otherwise % 'on'
   % Freezes:
 
   % Gets colorbar axes properties:
   cb_prop  = get(h);
   
   % Gets colorbar image handle. Fixed BUG, Sep 2009
   hi = findobj(h,'Type','image');
    
   % Gets image data and transform it in a RGB:
   CData = get(hi,'CData'); 
   if size(CData,3)~=1
    % It's already frozen:
    continue
   end
  
   % Gets image tag:
   Tag = get(hi,'Tag');
  
   % Deletes previous colorbar preserving peer axes position:
   oldunits = get(pax,'Units');
              set(pax,'Units','Normalized')
   Position = get(pax,'Position');
   delete(h)
   cbfreezedata.Position = get(pax,'Position');
              set(pax,'Position',Position)
              set(pax,'Units',oldunits)
  
   % Generates new colorbar axes:
   % NOTE: this is needed because each time COLORMAP or CAXIS is used,
   %       MATLAB generates a new COLORBAR! This eliminates that behaviour
   %       and is the central point on this function.
   h = axes(...
    'Parent'  ,cb_prop.Parent,...
    'Units'   ,'Normalized',...
    'Position',cb_prop.Position...
   );
  
   % Save location for future call:
   cbfreezedata.Location = cb_prop.Location;
  
   % Move ticks because IMAGE draws centered pixels:
   %XLim = cb_prop.XLim;
   %YLim = cb_prop.YLim;
   %if     isempty(cb_prop.XTick)
    % Vertical:
   % X = XLim(1) + diff(XLim)/2;
   % Y = YLim    + diff(YLim)/(2*length(CData))*[+1 -1];
   %else % isempty(YTick)
    % Horizontal:
   % Y = YLim(1) + diff(YLim)/2;
   % X = XLim    + diff(XLim)/(2*length(CData))*[+1 -1];
   %end
  
   % Draws a new RGB image:
   image(X,Y,ind2rgb(CData,colormap),...
    'Parent'            ,h,...
    'HitTest'           ,'off',...
    'Interruptible'     ,'off',...
    'SelectionHighlight','off',...
    'Tag'               ,Tag...
   )  

   % Removes all   '...Mode'   properties:
   cb_fields = fieldnames(cb_prop);
   indmode   = strfind(cb_fields,'Mode');
   for k=1:length(indmode)
    if ~isempty(indmode{k})
     cb_prop = rmfield(cb_prop,cb_fields{k});
    end
   end
   
   % Removes special COLORBARs properties:
   cb_prop = rmfield(cb_prop,{...
    'CurrentPoint','TightInset','BeingDeleted','Type',...       % read-only
    'Title','XLabel','YLabel','ZLabel','Parent','Children',...  % handles
    'UIContextMenu','Location',...                              % colorbars
    'ButtonDownFcn','DeleteFcn',...                             % callbacks
    'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle',...
    'PlotBoxAspectRatio','DataAspectRatio','Position',... 
    'XLim','YLim','ZLim'});
   
   % And now, set new axes properties almost equal to the unfrozen
   % colorbar:
   set(h,cb_prop)

   % CallBack features:
   set(h,...
    'ActivePositionProperty','position',...
    'ButtonDownFcn'         ,@cbfreeze,...  % mhh...
    'DeleteFcn'             ,@cbfreeze)     % again
   set(pax,'DeleteFcn'      ,@cbfreeze)     % and again!  
  
   % Do not zoom or pan or rotate:
   setAllowAxesZoom  (zoom    ,h,false)
   setAllowAxesPan   (pan     ,h,false)
   setAllowAxesRotate(rotate3d,h,false)
   
   % Updates data:
   CBH(icb) = h;   

   % Saves data for future undo:
   cbfreezedata.pax       = pax;
   setappdata(  h,cbappname,cbfreezedata);
   setappdata(pax,axappname,h);
   
 end % switch functionality   

end  % MAIN loop


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Output?:
if ~nargout
 clear CBH
else
 CBH(~ishandle(CBH)) = [];
end

% Returns current axes:
if ishandle(cax) 
 set(fig,'CurrentAxes',cax)
end

end % [EOF]   cbfreeze.m



function CBH = cbhandle(varargin)
%CBHANDLE   Handle of current colorbar axes.
%
%   SYNTAX:
%     CBH = cbhandle;
%     CBH = cbhandle(H);
%
%   INPUT:
%     H - Handles axes, figures or uipanels to look for colorbars.
%         DEFAULT: gca (current axes)
%
%   OUTPUT:
%     CBH - Color bar handle(s).
%
%   DESCRIPTION:
%     By default, color bars are hidden objects. This function searches for
%     them by its 'axes' type and 'Colorbar' tag.
%    
%   SEE ALSO:
%     COLORBAR
%     and
%     CBUNITS, CBLABEL, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbhandle.m
%   VERSION: 1.1 (Aug 20, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.1      Fixed bug with colorbar handle input. (Aug 20, 2009)

%   DISCLAIMER:
%   cbhandle.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
axappname = 'FrozenColorbar'; % Peer axes application data with frozen
                              % colorbar handle.

% Sets default:
H = get(get(0,'CurrentFigure'),'CurrentAxes');

if nargin && ~isempty(varargin{1}) && all(ishandle(varargin{1}))
 H = varargin{1};
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Looks for CBH:
CBH = [];
% set(0,'ShowHiddenHandles','on')
for k = 1:length(H)
 switch get(H(k),'type')
  case {'figure','uipanel'}
   % Parents axes?:
   CBH = [CBH; ...
    findobj(H(k),'-depth',1,'Tag','Colorbar','-and','Type','axes')];
  case 'axes'
   % Peer axes?:
   hin  = double(getappdata(H(k),'LegendColorbarInnerList'));
   hout = double(getappdata(H(k),'LegendColorbarOuterList'));
   if     (~isempty(hin)  && ishandle(hin))
    CBH = [CBH; hin];
   elseif (~isempty(hout) && ishandle(hout))
    CBH = [CBH; hout];
   elseif isappdata(H(k),axappname)
    % Peer from frozen axes?:
    CBH = [CBH; double(getappdata(H(k),axappname))];
   elseif strcmp(get(H(k),'Tag'),'Colorbar') % Fixed BUG Aug 2009
    % Colorbar axes?
    CBH = [CBH; H(k)];
   end
  otherwise
   % continue
 end
end
% set(0,'ShowHiddenHandles','off')

end % [EOF]   cbhandle.m

