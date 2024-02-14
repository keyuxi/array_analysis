

function AnalyseImage(seqDataFilename, imageFilename, dataScaling, filterSubsets, registrationOffsetMapFilename, maskImageFilename,CPfluorFilename)
    %### DEBUG ###    
    if(GlobalVars.debug)
        global debugFigMan;
        debugFigMan = debugFigureManager();
    end
    %### END DEBUG ###

    %read in image
    disp('Loading image...')
    rawImage = double(imread(imageFilename));
       
    %dissect the image filename for naming various output files
    [imagePath, imageName, imageExt] = fileparts(imageFilename);
    if(~isempty(imagePath))
        imagePath = [imagePath StringFun.pathSlash()];
    end
    
    %get output filename
    if(~exist('CPfluorFilename','var'))
        outputPath = imagePath;
        outputName = imageName;
        outputExt = '.CPfluor';
        CPfluorFilename = [outputPath outputName outputExt]; %create output filename for fit, quantified fluorescence values with a unique file extension
    else
        [outputPath, outputName, outputExt] = fileparts(CPfluorFilename);
        if(~isempty(outputPath))
            outputPath = [outputPath StringFun.pathSlash()];
        end
        if(strcmpi(outputExt,'CPfluor'))
            error('AnalyseImage:AnalyseImage','Invalid output filename "%s". The file extension of the output file must be ".CPfluor"',CPfluorFilename);
        end
    end
    
    %load sequence data
    disp('Loading and scaling sequence data...')
    cDat = FileFun.loadCPseq(seqDataFilename, dataScaling);

    %load the precomputed registration offset map from file
    registrationOffsetMap = RegistrationFun.loadRegistrationOffsetMap(registrationOffsetMapFilename);

    %for MiSeq/ImagingStation images, load the field of view mask image to eliminate data points outside the circular field of illumination
    if(strcmp(dataScaling,'MiSeq_to_ImagingStation'))
        maskImage = RegistrationFun.loadFOVmaskImage(maskImageFilename);
        %### DEBUG ###
        if(GlobalVars.debug)
            showDebug = debugFigMan.debugStop('DEBUG:AnalyseImage: Would you like to see the mask image?');
            if(showDebug)
                msg = sprintf('DEBUG:AnalyseImage: Mask image from file "%s"', maskImageFilename);
                PlotFun.plotClusterImage(maskImage, msg);
            end
            drawnow;
        end
        %### END DEBUG ###
    else
        maskImage = [];
    end

    %mark the subset of clusters that will be used for registration
    cDat.currSubset = RegistrationFun.identifyFilterSubsetsForRegistration(cDat, filterSubsets);
    
    %do the initial global registration
    bgImage = ImageFun.getBackgroundFromClusterImage(rawImage);% remove background from the image with morphological opening (this greatly improves cross correlation performance)
    bgSubtractedImage = rawImage - bgImage;
    disp('Initial global registration...');
    %### DEBUG ###
    if(GlobalVars.debug)
        showDebug = debugFigMan.debugStop('DEBUG:AnalyseImage: Would you like to get global offset by cross correlating to another image?');
        if(showDebug)
            
            % load reference image to cross correlate to
            [refImageFile, refImagePath] = uigetfile('*.tif');
            refImage = double(imread(fullfile(refImagePath, refImageFile)));
            bgRefImage = ImageFun.getBackgroundFromClusterImage(refImage);% remove background from the image with morphological opening (this greatly improves cross correlation performance)
            bgSubtractedRefImage = refImage - bgRefImage;
            
            % find registration offset of this image
            [cDat,globalOffset_y,globalOffset_x] = RegistrationFun.RegisterGlobal(bgSubtractedRefImage, cDat);
            
            % find offset between this image and the actual image to analyse
            [offsetY, offsetX, registrationSuccess, correlationPeakHeight] = RegistrationFun.RegisterTileByImage(bgSubtractedImage, bgSubtractedRefImage);
            if registrationSuccess;
                cDat.adjY = cDat.adjY - offsetY;
                cDat.adjX = cDat.adjX - offsetX;
            else
                sprintf('DEBUG:AnalyseImage: Image cross correlation was unsuccessful');
            end
        else
            [cDat,globalOffset_y,globalOffset_x] = RegistrationFun.RegisterGlobal(bgSubtractedImage, cDat);    
        end
    else
        [cDat,globalOffset_y,globalOffset_x] = RegistrationFun.RegisterGlobal(bgSubtractedImage, cDat);
    end
    
    %### DEBUG ###
    if(GlobalVars.debug)
        disp(sprintf('\n      The global y offset is "%f"', globalOffset_y))
        disp(sprintf('      The global x offset is "%f"', globalOffset_x))
        
        showDebug = debugFigMan.debugStop('DEBUG:AnalyseImage: Would you like to see the initial global registration of sequence data to the image?');
        if(showDebug)
            PlotFun.plotGlobalRegistrationByFilter(rawImage, cDat, [], [], 'initial global registration');
        end
        drawnow;
    end
    %### END DEBUG ###  
    
    %apply the pre-computed registration offset map and eliminate clusters outside of the illuminated image region
    disp('Applying the fine registration offset map to the globally registered points...');
    [cDat,localOffsets_y,localOffsets_x] = RegistrationFun.ApplyRegistrationOffset(rawImage,cDat,registrationOffsetMap,size(rawImage),maskImage);

    %### DEBUG ###
    if(GlobalVars.debug)
        showDebug = debugFigMan.debugStop('DEBUG: AnalyseImage: Would you like to see the final registration of sequence data to the image?');
        if(showDebug)
            PlotFun.plotGlobalRegistrationByFilter(rawImage, cDat, [], [], 'final registration (excised)');
        end
        
        %crop image to work on a smaller subtile for debugging
        response = DisplayFun.txtmenu('DEBUG:AnalyseImage:AnalyseImage: crop image to a smaller central subtile for debugging the quanitification?','crop to 300x300 px','use full image');
        if(response==0) %crop to 300x300
            [rawImage, offset] = ImageFun.getCentralSubtileImage(rawImage, 300, 300);
            cDat.adjY = cDat.adjY - offset(1); %offset all cluster positions to account for cropping
            cDat.adjX = cDat.adjX - offset(2);
        end
        drawnow;
    end
    %### END DEBUG ###
    
    % Clean image (i.e. remove bright blobs caused by fluorescent protein aggregation on the chip) %%%
    disp('Cleaning image...')
    % cleanimage looks for pixel values above maxVal, then removes the bright region immediately around it until it falls below cutOff
    % maxVal = 4000; % Could use 2*allclusterimage
    % cutOff = 1000;
    %[I, zeroPix] = CleanImage_v8(rawImage,maxVal,cutOff);
    I = rawImage;
    zeroPix = [];   


    %quantify image
    disp('Quantifying image...')
    [cDat,fitBgImage] = QuantifyImage(I,cDat,zeroPix);

    %save quantified fluorescence values to a text file HKWS: moved this up
    disp(['Saving quantified flourescence data to file: ' CPfluorFilename]);
    FileFun.saveCPfluorFile(cDat,CPfluorFilename);

    %generate and save QC image
    if(GlobalVars.makeQCfigs || GlobalVars.debug)
        qcPath = [outputPath 'QCFits' StringFun.pathSlash()];
        if(~exist(qcPath,'dir'))
            mkdir(qcPath);
        end
        qcFilename = [qcPath imageName '.FitQC.eps']; %create output filename for fit QC figure        
        disp(['Saving fit QC figure to file: ' qcFilename]);
        PlotFun.saveThreePanelFitWithoutGui(I,cDat,fitBgImage,[],[],qcFilename,globalOffset_y,globalOffset_x);
    end
    
    %### DEBUG ###
    if(GlobalVars.debug)
        showDebug = debugFigMan.debugStop('DEBUG: AnalyseImage: Would you like to see the final 3-panel image fit?');
        if(showDebug)
            %plot QC image
            PlotFun.plotThreePanelFit(I, cDat, fitBgImage);
        end
    end
    

    %done
    disp('Image analysis complete...');
