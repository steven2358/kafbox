% This is the template used for kernel adaptive filtering algorithms in
% the kernel adaptive filtering toolbox. Delete this line.
%
% The name of the algorithm goes here.
%
% A reference to the original publication of the algorithm goes here.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef kafbox_template
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        param1 = 1;
        param2 = 2;
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'protected', SetAccess = 'private') % variables
        dict = []; % dictionary
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = kafbox_template(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if strmatch(fn,fieldnames(kaf),'exact'),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kaf.alpha = 0;
            else
                
                % main algorithm training goes here
                
                % example of a helper function
                kaf = helper1(kaf,x,y);
            end
            
        end
        
    end
    
    methods (Access = 'private') % helper functions go here
        
        function kaf = helper1(kaf,x,y)
            % operations
        end
        
    end
end
