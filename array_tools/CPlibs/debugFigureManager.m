classdef debugFigureManager < handle
    
    properties
        locallyDismissed = containers.Map(); %keeps track of debug stops that have been locally dismissed and will not be shown again during local execution, but will come back the next time around
        globallyDismissed = containers.Map(); %keeps track of debug stops that have been globally dismissed and will not be shown anymore ever
    end % END properties
    
    methods % public methods

        %constructor
        function this = debugFigureManager()
            remove(this.locallyDismissed,keys(this.locallyDismissed)); %clear lists every time it is instantiated
            remove(this.globallyDismissed,keys(this.globallyDismissed));
        end
        
        function showDebug = debugStop(this, message)
            showDebug = false;
            if(~isKey(this.globallyDismissed,message) && ~isKey(this.locallyDismissed,message))
                response = DisplayFun.txtmenu(message,'yes','no more for now','no more of this ever');
                if(response==0)
                    showDebug = true;
                elseif(response==1)
                    this.locallyDismissed(message) = length(dbstack);
                elseif(response==2)
                    this.globallyDismissed(message) = true;
                end
            end
        end %end localDebug
        
        function clearLocal(this)
            currKeys = keys(this.locallyDismissed);
            currStackLevel = length(dbstack);
            for i = 1:length(currKeys)
                if(this.locallyDismissed(currKeys{i}) >= currStackLevel)
                    remove(this.locallyDismissed,currKeys{i});
                end
            end
        end
        
        function delete(this)
            remove(this.locallyDismissed,keys(this.locallyDismissed)); %clear lists every time it is deleted
            remove(this.globallyDismissed,keys(this.globallyDismissed));
        end
        
    end % END public methods
    
end

