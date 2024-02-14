% Container class that groups together the arrays that store all the information about the clusters in (typically) one tile.

classdef ClusterData < handle
    
    properties
        % these should all be equal-length arrays, where each entry of each
        % array represents the data for one cluster.
        IDs % unique ID for each cluster, taken from a combination of the machine ID, chip ID and cluster x,y position
        filterMatches; %colon-delimited string of filter names (alignment-based classifiers) that this clusters' sequence(s) match
        origX % original cluster x,y positions from data file. These should not be changed.
        origY
        
        numClusters; %total number of clusters
        
        excised % instead of throwing out clusters that we do not want to process, we mark them in this array (if entry == 1, then this cluster is excised and should not be processed)
        alreadyMoved; %temporary flag to keep track of clusters that have already been moved as clusters are moved in regions
        alreadyFit; %flag to keep track of clusters that have already been fit to image data 
        
        adjX % adjusted cluster x,y positions (adjusted to account for correct registration of the image)
        adjY
        
        %fit parameters
        fitAmplitudes;
        fitY;
        fitX;
        fitSigmas;
        fitBackground;

        ID_indexMap; %maps the unique IDs to the indicdes of the clusters
        
        filterSubsets; %a map of lists of the indices of each filtered subset, keyed by the filter names
        
        currSubset;  %logical array incicating which clusters are part of the current operation (e.g. have not excised and match certain filters of interest)
    end
    
    methods %public methods
        function this = ClusterData(init_IDs, init_filterMatches, filterMap, init_x, init_y) %constructor
            if(~isequal(size(init_IDs),size(init_filterMatches),size(init_x),size(init_y)))
                error('ClusterData:ClusterData:unequalDataLenghts','All cluster data must be in identically indexed arrays of the same length');
            end
            
            this.IDs = init_IDs;
            this.filterMatches = init_filterMatches;
            this.origX = init_x;
            this.origY = init_y;
            
            this.numClusters = length(init_x);
            
            this.excised = false(this.numClusters,1);
            this.alreadyMoved = false(this.numClusters,1);
            this.alreadyFit = false(this.numClusters,1);
            
            this.adjX = this.origX;
            this.adjY = this.origY;
            
            this.fitAmplitudes = zeros(this.numClusters,1);
            this.fitY = zeros(this.numClusters,1);
            this.fitX = zeros(this.numClusters,1);
            this.fitSigmas = zeros(this.numClusters,1);
            this.fitBackground = zeros(this.numClusters,1);

            this.ID_indexMap = containers.Map(this.IDs,1:this.numClusters); %a map of cluster IDs to indices 
            
            %store a map of logical arrays indicating which elements belong to each filtered subset
            this.filterSubsets = containers.Map();
            filterNameArray = keys(filterMap);
            matchesAnyCluster = false(this.numClusters,1);
            for i = 1:size(filterNameArray,2)
                matchesThisFilter = ~cellfun(@isempty,strfind(this.filterMatches,filterNameArray{i})); %logic indicates if the cluster matches the filter
                this.filterSubsets(filterNameArray{i}) = matchesThisFilter; %add to the map
                matchesAnyCluster = matchesAnyCluster | matchesThisFilter;
            end
            this.filterSubsets('__Matches_no_filters__') = ~matchesAnyCluster; %a final category for clusters that don't match any filter
            
            this.currSubset = false(this.numClusters,1); %logical array incicating which clusters are part of the current operation (e.g. have not excised and match certain filters of interest)
        end

        function index = getIndex(this, ID)
            index = this.ID_indexMap(ID);
        end
        
        function tf = isValidFilterName(this, filterName)
            tf = isKey(this.filterSubsets, filterName);
        end
        
        function subset = getFilteredSubset(this, filterName) %returns a logical array
            subset = [];
            if isKey(this.filterSubsets, filterName)
                subset = this.filterSubsets(filterName);
            end               
        end
        
        function x = currX(this)
            x = this.adjX(this.currSubset);
        end
        
        function y = currY(this)
            y = this.adjY(this.currSubset);
        end
        
        function resetAlreadyMoved(this)
            this.alreadyMoved = false(this.numClusters,1);
        end
        
        function updateExcised(this, toBeExcised) %pass in a logical array of the clusters to be excised
            this.excised = this.excised | toBeExcised;
            this.currSubset = this.currSubset & ~toBeExcised;
        end
    end % END public methods
end
