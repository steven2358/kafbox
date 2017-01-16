% Fixed-budget kernel least mean squares (FB-KLMS) algorithm.
%
% D. Rzepka, "Fixed-budget kernel least mean squares," 2012 IEEE 17th
% Conference on Emerging Technologies & Factory Automation (ETFA), Krakow,
% Poland, Sept. 2012, http://dx.doi.org/10.1109/ETFA.2012.6489767
%
% Remark: code contributed by Dominik Rzepka
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef fbklms < handle
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        nu = .05; % growth criterion threshold
        M = 500; % dictionary size
        eta = .5; % learning rate
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'private', SetAccess = 'private') % variables
        dict = []; % dictionary
        diagkdict = []; % diagonal of kernel matrix for dictionary
        alpha = []; % expansion coefficients
    end
    
    methods        
        function kaf = fbklms(parameters) % constructor
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
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                kaf.diagkdict(1) = k;
                kaf.alpha = kaf.eta*y*k/(k'*k);
            else
                
                kt = kernel([kaf.dict;x],x,kaf.kerneltype,kaf.kernelpar);
                k = kt(1:end-1);
                y_est = k'*kaf.alpha;
                e = y - y_est;
                
                kaf.alpha = kaf.alpha + kaf.eta*e*k/(k'*k);
                
                % growth criterion
                kx = kt(end);
                dependency = kx - k./kaf.diagkdict;
                if min(dependency) >= kaf.nu, % expand dictionary
                    kaf.dict = [kaf.dict; x];
                    kaf.diagkdict = [kaf.diagkdict; kx];
                    kaf.alpha = [kaf.alpha; kaf.eta*e/(k'*k)];
                    
                    if length(kaf.alpha) > kaf.M, % prune dictionary
                        [~, id] = min(abs(kaf.alpha));
                        kaf.dict(id,:) = [];
                        kaf.diagkdict(id) = [];
                        kaf.alpha(id) = [];
                    end
                end
            end
        end
    end
end
