% Kernel Least-Mean-Square algorithm
% Author: Steven Van Vaerenbergh, 2013
% Reference: http://dx.doi.org/10.1109/TSP.2007.907881
% Comment: implementation includes a maximum dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef klms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu = 0.5; % learning rate
        M = 1000; % maximum dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mem = []; % memory
        alpha = []; % expansion coefficients
    end
    
    methods
        
        function kaf = klms(parameters) % constructor
            if (nargin > 0)
                kaf.mu = parameters.mu;
                kaf.M = parameters.M;
                kaf.kerneltype = parameters.kerneltype;
                kaf.kernelpar = parameters.kernelpar;
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.mem,1)>0
                k = kernel(kaf.mem,x,kaf.kerneltype,kaf.kernelpar);
                y_est = k'*kaf.alpha;
            else
                y_est = 0;
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            y_est = kaf.evaluate(x);
            err = y - y_est;
            
            if (size(kaf.mem,1)<kaf.M),
                kaf.alpha = [kaf.alpha; kaf.mu*err]; % grow
                kaf.mem = [kaf.mem; x]; % grow
            end
        end
        
    end    
end
