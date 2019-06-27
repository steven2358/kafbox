% Kernel Recursive Least Squares algorithm with Approximate Linear
% Dependency criterion
%
% Y. Engel, S. Mannor, and R. Meir, "The kernel recursive least-squares
% algorithm," IEEE Transactions on Signal Processing, vol. 52, no. 8, pp.
% 2275-2285, Aug. 2004, http://dx.doi.org/10.1109/TSP.2004.830985
%
% Remark: implementation includes a maximum dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef krls < kernel_adaptive_filter
    
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
    end
    
    methods
        
        function kaf = krls(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if ismember(fn,fieldnames(kaf)),
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
        
        function train(kaf,x,y) % train the algorithm
            k = kernel([kaf.dict; x],x,kaf.kerneltype,kaf.kernelpar);
            kt = k(1:end-1);
            ktt = k(end);
            if numel(kt)==0, % initialize
                kaf.Kinv = 1/ktt;
                kaf.alpha = y/ktt;
                kaf.P = 1;
                kaf.dict = x;
            else
                at = kaf.Kinv*kt; % coefficients of best linear combination
                delta = ktt - kt'*at; % residual of linear approximation
                
                if (delta>kaf.nu && size(kaf.dict,1)<kaf.M), % not ALD
                    kaf.dict = [kaf.dict; x]; % add base to dictionary
                    kaf.Kinv = 1/delta*... % update inverse kernel matrix
                        [delta*kaf.Kinv + at*at', -at; -at', 1];
                    Z = zeros(size(kaf.P,1),1);
                    kaf.P = [kaf.P Z; Z' 1]; % extend projection matrix
                    ode = 1/delta*(y-kt'*kaf.alpha);
                    kaf.alpha = [kaf.alpha - at*ode; ode]; % full update
                else % perform reduced update of coefficients
                    q = kaf.P*at/(1+at'*kaf.P*at);
                    kaf.P = kaf.P - q*(at'*kaf.P); % update proj. matrix
                    kaf.alpha = kaf.alpha + kaf.Kinv*q*(y-kt'*kaf.alpha);
                end
            end
        end
    end
end
