% This is the template used for kernel adaptive filtering algorithms in
% the kernel adaptive filtering toolbox.
%
% The name of the algorithm goes here.
%
% A reference to the original publication of the algorithm goes here.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef kaf_template % replace this with the algorithm acryonym
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters with their default values
        param1 = 1;
        param2 = 2;
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % internal variables
        dict = []; % dictionary
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = kaf_template(parameters) % constructor, replace this with the algorithm acryonym
            if (nargin > 0) % replace with parameters
                kaf.param1 = parameters.param1;
                kaf.param2 = parameters.param2;
                kaf.kerneltype = parameters.kerneltype;
                kaf.kernelpar = parameters.kernelpar;
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if size(kaf.dict,2)==0 % initialize if necessary
                kaf.dict = x;
                kaf.alpha = 0;
            else
                
                % main algorithm training goes here
                
                kaf = helper1(kaf,x,y); % example of a helper function to keep code modular
            end
            
        end
        
    end
    
    methods (Access = 'private') % helper functions go here
        
        function kaf = helper1(kaf,x,y)
            % operations
        end
        
    end
end
