classdef StringFun
    methods (Static)
        
        function s = pathSlash() %returns the appropriate slash for OS-specific path names
            if(ispc()) %if is a windows os
                s = '\';
            else %unix or mac
                s = '/';
            end
        end
        
         function [cellArrayOfTokens, numTokens] = tokenizeString(inString, delimiter)
            c = textscan(inString,'%s','delimiter',delimiter);
            c = c{1};
            numTokens = length(c);
            for i = 1:numTokens
                cellArrayOfTokens{i} = c{i};
            end
        end
        
    end % END static methods
end