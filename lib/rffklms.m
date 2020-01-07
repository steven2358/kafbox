% Random Fourier Feature Kernel Least Mean Square algorithm.
%
% Abhishek Singh, Narendra Ahuja and Pierre Moulin, "Online learning with
% kernels: Overcoming the growing sum problem," 2012 IEEE International
% Workshop on Machine Learning for Signal Processing (MLSP), Sept. 2012.
% http://dx.doi.org/10.1109/MLSP.2012.6349811
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef rffklms < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        seed = 1;
        mu = .9; % step size
        D = 1000; % RFF dimension
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        omega = []; %
        Omega = []; % weight vector
        b = 0; %
    end
    
    methods        
        function kaf = rffklms(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)'
                    if ismember(fn,fieldnames(kaf))
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if numel(kaf.omega)
                Psi = cos(kaf.omega*x' + repmat(kaf.b,1,size(x,1)));
                y_est = Psi'*kaf.Omega/kaf.D;
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function train(kaf,x,y) % train the algorithm
            if ~numel(kaf.omega)
                rng('default');
                rng(kaf.seed);
                kaf.omega = 1/kaf.kernelpar*randn(kaf.D,size(x,2));
                kaf.b = 2*pi*rand(kaf.D,1);
                kaf.Omega = zeros(kaf.D,1);
            end
            
            Psi = cos(kaf.omega*x' + kaf.b);
            y_est = kaf.Omega'*Psi/kaf.D;
            err = y - y_est;
            kaf.Omega = kaf.Omega + kaf.mu*err*Psi;
        end        
    end
end
