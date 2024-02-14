% Global vars for image analysis/data processing pipeline

classdef GlobalVars
     properties (Constant)
         
         debug = false; %debug mode outputs debugging information and figures
         
         %generate QC images for the fits and save as .eps files?
         makeQCfigs = true;
         
         %scaling factors for sequencing data generated on one machine then
         %imaged on another
         scalingFactor = containers.Map({'MiSeq_to_ImagingStation', 'MiSeq_to_TIRFStation1', 'GA_to_GA', 'none'} , {10.9240, 11.13, 10, 1});         
                 
         %when making a synthetic image of gaussians, each individual
         %gaussian is computed inside of a local window defined by this param.
         %The tails of a normal distribution go to infinity, but for
         %computational efficiency, we compute each gaussian in a local
         %window, the size of which is chosen so that no signal is lost 
         %that is higher than the "AmplitudeCutoffForLocalWindow" factor
         %below the peak amplitude
         AmplitudeCutoffForLocalWindow = 1e15;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Registration Parameters %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %pixel intensity value for data points on a synthetic image used in cross-correlation
         synthIntensityValue = 10;
         
         %if the crosscorrelation peak height is at least this value, registration (global or of a subtile) is considered to have succeeded 
         correlationSuccessPeakHeight = 0.025;
         
         %default width and height of the central subtile that is used for global registration
         globalRegistrationSubtileHeight = 500;
         globalRegistrationSubtileWidth = 500;
         
         %maxiuum offset from the previous position that is allowed in global registration
         %maxOffsetGlobalRegistration = [200 200];
         maxOffsetGlobalRegistration = [1000 500];

         
         %number of tiles used at each levels of hierarchical registration
         numTilesCoarseHierarchicalRegistration = [4 4];
         numTilesFineHierarchicalRegistration = [16 16];
         
         %maxiuum offset from the previous position that is allowed at each level of hierarchical registration
         maxOffsetCoarseHierarchicalRegistration = [70 70];
         maxOffsetFineHierarchicalRegistration = [20 20];
         
         %bound around the initial putative pixel-resolution correlation peak that in wich we constrain the search for the sub-pixel peak of the gaussian
         corrPeak_peakBound = 11.5;
         
         %data window around the initial putative correlation peak that fit with a gaussian-fit correlation peak (in pixels)
         corrPeak_window = 7;
         
         %range of sigma values in fiting the correlation peak to a single 2D gaussian [lb initialGuess ub]
         corrPeak_sigmaRange = [0.25 1.42 13.00];
         
         %tolerances and limits for the least-squares fiting of the correlation peak to a single 2D gaussian 
         corrPeak_TolX = 1e-4;
         corrPeak_TolFun = 1e-11;
         corrPeak_MaxFunEvals = 100;
         
         %if a subtile contains fewer than this number of clusters, we
         %consider it unregisterable
         minNumClustersForRegistration = 15;
         
         %bounding box around the correlation peak used to calculate background
         corrPeak_background_window = 10;

         %checkTileRegistration - parameters for finding registration by registering multiple subtiles at lower thresholds to the same frame of registration
         matchingCorrelationTolerance = 15; %subtile registrations must be within this number of pixels in both axes to be considered in the same frame
         numTileToMatch = 2; %number of subtiles that must match at a lower threshold factor in the same frame of registration to be considered a match
         nTilesRegistrationCheck = [3 3]; %subtile matrix for checkRegistration
         cropBorder = 241; %border to crop away around each subtile (in pixels)
         corrSuccessPeakHeightRegistrationCheck = 0.1; %different threshold for correlation success when using checkRegistration
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Auto FOV Mask Generation Parameters %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %when auto-generating a field of view mask, the background illumination
         %pattern is smoothed with a disk filter of this radius
         smoothingFilterRadius = 20;
         
         %when auto-generating a field of view mask, regions with
         %background illumination that are not at least this fraction of
         %the mean intensity will be considered unilluminated and eliminated
         maskIntensityCutoff = 0.60;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Gaussian Image Fitting Parameters %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %the image will be fit to multiple gaussians in small overlapping subtiles
         
         %the size (in pixels) of the inner window
         innerWindowSize = 18;
         
         %the size (in pixels) of the overlap between subtiles
         overlap = 4;
         
         %the amount that the individual cluster positions are allowed to move from the registred position during fitting
         positionTolerance = 0.3;
         
         %range of sigma values in the fit optimization [lb initialGuess ub]
         multiGauss_sigmaRange = [0.95 1.29 1.55];
         
         %tolerances and limits for the least-squares optimization
         multiGauss_TolX = 0.05;
         multiGauss_TolFun = 1e-6;
         multiGauss_MaxFunEvals = 200;
         
         %%%%%%%%%%%%%%%
         %%% File IO %%%
         %%%%%%%%%%%%%%%
         
         chunkSize = 1000; %large text files are read in/out in chunks of this many lines at a time.  MATLAB's slow file access times that make per-line file access very slow

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Sequence Filtering Parameters %%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %statistical parameters that were fit to the distribution of scores of randomly-generated sequences 20 - 40bp long   
         EVD_k = 0.0097; %constant emperically fit to an extreme value distribution
         EVD_lambda = 0.5735; %constant, empirically fit to an extreme value distribution

         pThreshold = 1e-4; % threshold for NW alignment to be accepted as "sufficiently close" that it passes the filter

         %%%%%%%%%%%%%%%%%%%%
         %%% Random Stuff %%%
         %%%%%%%%%%%%%%%%%%%%
         
         goldenMean = 1.618; % visually appealing ratio for setting intensity thresholds for image display
                  
     end % END static properties
end % END GlobalVars class
