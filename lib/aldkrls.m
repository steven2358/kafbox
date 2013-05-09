% Kernel Recursive Least Squares algorithm with Approximate Linear
% Dependency
% Author: Steven Van Vaerenbergh, 2013
% Reference: http://dx.doi.org/10.1109/TSP.2004.830985
% Comment: implementation includes a maximum dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef aldkrls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        nu = 1E-4; % ALD threshold
        M = 1000; % maximum dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        alpha = []; % expansion coefficients
        P = []; % inverse A'*A
        Kinv = []; % inverse kernel matrix
        grow = false; % flag
    end
    
    methods
        
        function kaf = aldkrls(parameters) % constructor
            if (nargin > 0)
                kaf.nu = parameters.nu;
                kaf.M = parameters.M;
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
            k = kernel([kaf.dict; x],x,kaf.kerneltype,kaf.kernelpar);
            kt = k(1:end-1);
            ktt = k(end);
            if numel(kt)==0 % initialize
                kaf.Kinv = 1/ktt;
                kaf.alpha = y/ktt;
                kaf.P = 1;
                kaf.dict = x;
            else
                at = kaf.Kinv*kt;
                delta = ktt - kt'*at;
                
                if (delta>kaf.nu && size(kaf.dict,1)<kaf.M) % expand dictionary
                    kaf.grow = true;
                    kaf.dict = [kaf.dict; x];
                    kaf.Kinv = 1/delta*[delta*kaf.Kinv + at*at', -at; -at', 1];
                    Z = zeros(size(kaf.P,1),1);
                    kaf.P = [kaf.P Z; Z' 1];
                    ode = 1/delta*(y-kt'*kaf.alpha);
                    kaf.alpha = [kaf.alpha - at*ode; ode];
                else % only update alpha
                    kaf.grow = false;
                    q = kaf.P*at/(1+at'*kaf.P*at);
                    kaf.P = kaf.P - q*(at'*kaf.P);
                    kaf.alpha = kaf.alpha + kaf.Kinv*q*(y-kt'*kaf.alpha);
                end
            end
        end
        
    end
end
