classdef FitFun
    
    methods (Static)

        function [f,j] = quadraticSurface(params, yx, y0, x0)
            %assign params
            a1 = params(1);
            a2 = params(2);
            b1 = params(3);
            b2 = params(4);
            c = params(5);
            d = params(6);

            %retreive x,y vectors            
            y = yx(:,1);
            x = yx(:,2);

            %function value
            f = a1*((y-y0).^2) + a2*((x-x0).^2) + b1*(y-y0) + b2*(x-x0) + c*(y-y0).*(x-x0) + d;

            %jacobian
            if(nargout > 1)
                j(:,1) = (y-y0).^2;
                j(:,2) = (x-x0).^2;
                j(:,3) = y-y0;
                j(:,4) = x-x0;
                j(:,5) = (y-y0).*(x-x0);
                j(:,6) = ones(size(y));
            end
        end
        
        function [f,j] = gaussian_2D(params,yx)
            %assign params
            A = params(1);
            y0 = params(2);
            x0 = params(3);
            sigma = params(4);
            bg = params(5);

            %retreive x,y vectors            
            y = yx(:,1);
            x = yx(:,2);

            %function value
            expBase = exp(-1*(((x0-x).^2) + ((y0-y).^2)) / (2*(sigma^2)));
            f = bg + A*expBase; %add background
            
            %jacobian
            if(nargout > 1)
                j = ones(length(y),5); %initialize to ones
                j(:,1) = expBase;
                j(:,2) = -1*A*(y0-y).*(expBase/(sigma^2));  % y0
                j(:,3) = -1*A*(x0-x).*(expBase/(sigma^2));  % x0
                j(:,4) = A*(((x0-x).^2)+((y0-y).^2)).*(expBase/(sigma^3));  % sigma
                %j(:,5) = 1; %the derivative of bg is all ones, which is what we initialized to, so just we do not assign it
            end
        end

        %a function of N gaussians on a 2D grid (image)
        function [f,j] = multiGaussian_2D(params,yx)
            
            %init function return value
            f = 0;

            %retreive x,y vectors            
            y = yx(:,1);
            x = yx(:,2);
            
            if nargout > 1
                j = ones(length(y), length(params)); %if the jacobian is requested, initialize (to one)
            end
            
            %get global background parameter
            bg = params(end);
            
            % loop through clusters
            for i = 1:4:length(params)-4

                % get per-cluster gaussian parameter values
                A = params(i);
                y0 = params(i+1);
                x0 = params(i+2);
                sigma = params(i+3);
            
                %add function value to fit
                expBase = exp(-1*(((x0-x).^2) + ((y0-y).^2)) / (2*(sigma^2)));
                f = f + A*expBase;

                % jacobian
                if nargout > 1
                    j(:,i) = expBase;  % A
                    j(:,i+1) = -1*A*(y0-y).*(expBase/(sigma^2));  % y0
                    j(:,i+2) = -1*A*(x0-x).*(expBase/(sigma^2));  % x0
                    j(:,i+3) = A*(((x0-x).^2)+((y0-y).^2)).*(expBase/(sigma^3));  % sigma
                    %note that the derivative of the background is 1, which is what we initialized to, so we do not assign it
                end
            end
            
            % Add background
            f = f + bg;
        end
        
        
    end % END static methods
end