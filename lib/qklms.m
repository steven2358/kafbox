% Quantized Kernel Least Mean Square algorithm
% Author: Steven Van Vaerenbergh, 2013
% Reference: http://dx.doi.org/10.1109/TNNLS.2011.2178446
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef qklms
    
    properties (GetAccess = 'public', SetAccess = 'private')
        eta = .5; % learning rate
        epsu = 1;
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mem = []; % codebook
        alpha = []; % expansion coefficients
        grow = false; % flag
    end
    
    methods
        
        function kaf = qklms(parameters) % constructor
            if (nargin > 0)
                kaf.eta = parameters.eta;
                kaf.epsu = parameters.epsu;
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
          
            m = size(kaf.mem,1);
            if m==0
                d2 = kaf.epsu^2 + 1;
            else
                [d2,j] = min(sum((kaf.mem - repmat(x,m,1)).^2,2));
            end
            if d2 <= kaf.epsu^2, 
                kaf.grow = false;
                kaf.alpha(j) = kaf.alpha(j) + kaf.eta*err;
            else
                kaf.grow = true;
                kaf.mem = [kaf.mem; x]; % add to codebook
                kaf.alpha = [kaf.alpha; kaf.eta*err];
            end
        end
        
    end    
end
