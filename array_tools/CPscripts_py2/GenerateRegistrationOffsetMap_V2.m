%generates a registration offset map from sequencing data and an image of a
%given tile on a sequencing flow cell

%the input data used to calculate the map should be:
% 1.) an image (typically an all-cluster or densely-populated image) imaged in the same color and on the same machine as the data to be analyzed
% 2.) sequencing data sequenced on the same machine as the data to be analyzed

function GenerateRegistrationOffsetMap(seqDataFilename, imageFilename, dataScaling, filterSubsets, filterColormapFilename, registrationOffsetMapFilename, maskImageFilename)
    %### DEBUG ###
    if(GlobalVars.debug)
        global debugFigMan;
        if(isprop(debugFigMan,'globallyDismissed'))
            debugFigMan.delete();
        end
        debugFigMan = debugFigureManager();
    end
    %### END DEBUG ###
    
    %read in image
    disp(['Loading image from file ' imageFilename ':'])
    rawImage = double(imread(imageFilename));

    %### 20180711 Ben edit ###
    bgImage = ImageFun.getBackgroundFromClusterImage(rawImage);
    bgSubtractedImage = rawImage - bgImage;
    % ### use bgSubtractedImage for registration ###
   
    %dissect the image filename for naming various output files
    [imagePath, imageName, imageExt] = fileparts(imageFilename);
    if(~isempty(imagePath))
        imagePath = [imagePath StringFun.pathSlash()];
    end
    
    if(~exist('registrationOffsetMapFilename','var'))
        %default output filename for the registration offset map is the image name with a .roff extension
        registrationOffsetMapFilename = [imagePath imageName '.roff']; 
    end
    
    if(~exist('filterSubsets','var'))
        filterSubsets = []; %an empty filter subset name list means all clusters in the file will be used for registration
    end

    if(~exist('filterColormapFilename','var'))
        % User did not specify colormap for plots (to plot clusters matching different filters with different colors)
        % So we will specify a default colormap. The default depends on whether we are registering to a subset or not:
        % If we are registering to a subset, then we will define a colormap where each filter that is in the susbet is defined
        % a different color, and all filters not being used for registration are assigned to be blue.
        % If we are not registering to a subset, i.e. all clusters are being used for registration, we will allow not explicitly
        % define a colormap, and instead allow the plot function to assign colors (which it does by cycling through colors for each
        % filter).

        if ~isempty(filterSubsets)
	    filterColors = cell(length(filterSubsets),1);

	    availableColors = keys(PlotFun.colorMap);
	    % exclude 'blue', because we want to reserve that as a special color for all unused clusters
	    [tf, idx] = ismember('blue', availableColors);
	    if tf
                availableColors(idx) = []; % deletes 'blue' entry from cell array
            end

            for idx = 1:length(filterSubsets)
                filterColors{idx} = availableColors{ mod(idx, length(availableColors)) + 1 }; % assign each filter in the registration subset a color (reuse colors using MOD if there are more filters in the subset than there are available colors)
	    end
            filterColormap = containers.Map(filterSubsets, filterColors);
	    filterColormap('ELSE') = 'blue';
        else
            filterColormap = containers.Map(); % make empty Map; plot function will assign colors automatically
        end
    else
       filterColormap = FileFun.loadFilterColormap(filterColormapFilename);
    end
    
    %load sequence data
    disp('Loading and scaling sequence data...')
    cDat = FileFun.loadCPseq(seqDataFilename, dataScaling);
   
    %register image
    disp('Performing registration of the data to the image...')
    
    %Calculate the registration offset map from the current image with hierarchical registration
        
    if(strcmp(dataScaling,'MiSeq_to_ImagingStation'))
        %for MiSeq/ImagingStation images, automatically generate a mask image
        %to remove any data falling outside of the illuminated (circular) FOV
        
        if(~exist('maskImageFilename','var'))
            %default output filename for the FOV mask is the image name with a .mask.tif extension
            maskImageFilename = [imagePath imageName '.mask.tif']; %create output filename for registration offset map
        end  
        
        % create the field of view mask and save to to an image file 
        disp(['Generating FOV mask for MiSeq/ImagingStation images and saving to file: ' maskImageFilename]);
        maskImage = RegistrationFun.autoGenerateFOVmask(rawImage, maskImageFilename);
    else
        maskImage = [];
    end

    %mark the subset of clusters that will be used for registration
    cDat.currSubset = RegistrationFun.identifyFilterSubsetsForRegistration(cDat, filterSubsets);
    
    %generate registration offset map with hierarchical registration
    [registrationOffsetMap, globalOffset_y, globalOffset_x] = RegistrationFun.calculateRegistrationOffsetMap(bgSubtractedImage, cDat, maskImage, filterColormap);

    % save registration offset map to a text file
    disp(['Saving registration offset map to file: ' registrationOffsetMapFilename]);
    RegistrationFun.saveRegistrationOffsetMap(registrationOffsetMap, registrationOffsetMapFilename); %save the registration offset map in a text file
end
