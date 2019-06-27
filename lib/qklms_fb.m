% Fixed budget quantized kernel least-mean-square algorithm
%
% S. Zhao, B. Chen, P. Zhu, J. C. Principe, "Fixed budget quantized kernel
% least-mean-square algorithm", Signal Processing, Volume 93, Issue 9,
% September 2013, Pages 2759-2770,
% http://dx.doi.org/10.1016/j.sigpro.2013.02.012.
%
% Remark: significance calculation only implemented for Gaussian kernel.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef qklms_fb < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private')
        eta = .9; % learning rate
        epsu = .1; % quantization threshold
        beta = .95; % forgetting factor for influence
        M = 500; % dictionary size (K in publication)
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % codebook
        alpha = []; % expansion coefficients
        E = []; % significance vector (~ importance)
        lambda = []; % influence vector (~ how many samples are quantized
        % to each centre)
    end
    
    methods
        function kaf = qklms_fb(parameters) % constructor
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
            y_est = kaf.evaluate(x);
            err = y - y_est;
            
            m = size(kaf.dict,1);
            if m==0
                d2 = kaf.epsu^2 + 1; % force dictionary growth
            else
                [d2,j] = min(sum((kaf.dict - repmat(x,m,1)).^2,2));
            end
            
            c1 = pi*kaf.kernelpar^2/2; % pre-calculate
            
            if d2 <= kaf.epsu^2, % new basis under quantization threshold
                kaf.alpha(j) = kaf.alpha(j) + kaf.eta*err;
                
                % update significance for quantizing y to j-th centre (16)
                if m>1
                    inds = 1:m;
                    inds(j) = [];
                    kaf.E(inds) = kaf.beta*kaf.E(inds) + c1*abs(kaf.alpha(inds)).* ...
                        kernel(kaf.dict(inds,:),kaf.dict(j,:),kaf.kerneltype,kaf.kernelpar);
                end
                kaf.E(j) = abs(1 + kaf.eta*err/kaf.alpha(j)) * ...
                    kaf.beta*kaf.E(j) + ...
                    abs(kaf.alpha(j) + kaf.eta*err) * ...
                    c1*kernel(kaf.dict(j,:),kaf.dict(j,:),kaf.kerneltype,kaf.kernelpar);
                
                % update influence
                kaf.lambda = kaf.beta*kaf.lambda;
                kaf.lambda(j) = kaf.lambda(j) + 1;
                
            else % new basis not under quantization threshold
                if m < kaf.M % still room for extra centres
                    kaf.dict = [kaf.dict; x]; % add to codebook
                    kaf.alpha = [kaf.alpha; kaf.eta*err];
                    
                    % update significance for addition of m+1-th centre (15)
                    kaf.E(m+1,1) = 0;
                    kaf.E = kaf.beta*kaf.E + abs(kaf.alpha(m+1)) * ...
                        kernel(kaf.dict,kaf.dict(m+1,:),kaf.kerneltype,kaf.kernelpar);
                    
                    % update influence
                    kaf.lambda = kaf.beta*kaf.lambda;
                    % initial influence
                    kaf.lambda(m+1) = 1;
                    
                else % no room for extra centres
                    [~,L] = min(kaf.E); % centre with lowest significance
                    
                    % update significance for removing L-th centre (17)
                    kaf.E = kaf.E - kaf.lambda(L) * c1*abs(kaf.alpha) .* ...
                        kernel(kaf.dict,kaf.dict(L,:),kaf.kerneltype,kaf.kernelpar);
                    
                    % prune
                    kaf.dict(L,:) = [];
                    kaf.alpha(L) = [];
                    kaf.E(L) = [];
                    kaf.lambda(L) = [];
                    
                    % re-calculate error
                    y_est = kaf.evaluate(x);
                    err = y - y_est;
                    
                    % grow
                    kaf.dict = [kaf.dict; x];
                    kaf.alpha = [kaf.alpha; kaf.eta*err];
                    
                    % update significance for addition of m-th centre (15)
                    kaf.E(m,1) = 0;
                    kaf.E = kaf.beta*kaf.E + abs(kaf.alpha(m)) * ...
                        kernel(kaf.dict,kaf.dict(m,:),kaf.kerneltype,kaf.kernelpar);
                    
                    % update influence
                    kaf.lambda = kaf.beta*kaf.lambda;
                    % initial influence
                    kaf.lambda(m) = 1;
                end
            end
        end
    end
end
