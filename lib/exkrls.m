% Extended Kernel Recursive Least Squares algorithm
%
% W. Liu, I.M. Park. Y. Wang, and J.C. Principe, "Extended Kernel Recursive
% Least Squares Algorithm," IEEE Transactions on Signal Processing, vol.
% 57, no. 10, pp. 3801-3814, Oct. 2009,
% http://dx.doi.org/10.1109/TSP.2009.2022007
%
% Comment: implementation of the tracking model, includes a maximum
% dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef exkrls
    
    properties (GetAccess = 'public', SetAccess = 'private')
        alphaf = .999; % state forgetting factor, "alpha" in publication
        beta = .995; % data forgetting factor
        lambda = 1E-2; % regularization
        q = 1E-3; % trade-off between modeling variation and measurement disturbance
        M = 500; % maximum dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mem = []; % memory
        rho = [];
        Q = [];
        i = 0; % iteration number;
        alpha = []; % expansion coefficients, "a" in publication
        grow = false; % flag
    end
    
    methods
        
        function kaf = exkrls(parameters) % constructor
            allpars = {'alphaf','lambda','beta','q','kerneltype','kernelpar','M'};
            if (nargin > 0)
                for j=1:length(allpars),
                    p = allpars{j};
                    if isfield(parameters,p), kaf.(p) = parameters.(p); end
                end
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
            kaf.i = kaf.i + 1;
            k = kernel([kaf.mem; x],x,kaf.kerneltype,kaf.kernelpar);
            kt = k(1:end-1);
            ktt = k(end);
            kaf.grow = false;
            if numel(kt)==0 % initialize
                kaf.alpha = kaf.alphaf*y/(kaf.lambda*kaf.beta+ktt);
                kaf.rho = kaf.lambda*kaf.beta/(kaf.alphaf^2*kaf.beta + kaf.lambda*kaf.q);
                kaf.Q = kaf.alphaf^2/((kaf.beta*kaf.lambda+ktt)*(kaf.alphaf^2+kaf.beta*kaf.lambda*kaf.q));
                kaf.mem = x;
                kaf.grow = true;
            else
                if (size(kaf.mem,1)<kaf.M), % avoid infinite growth
                    z = kaf.Q*kt;
                    r = kaf.beta^kaf.i*kaf.rho + ktt - kt'*z;
                    err = y - kt'*kaf.alpha;
                    
                    kaf.alpha = kaf.alphaf*[kaf.alpha - z*err/r; err/r]; % grow
                    kaf.mem = [kaf.mem; x];
                    dummy = kaf.alphaf^2 + kaf.beta^kaf.i*kaf.q*kaf.rho;
                    kaf.rho = kaf.rho/dummy;
                    kaf.Q = kaf.alphaf^2/(r*dummy)*...
                        [kaf.Q*r + z*z', -z; -z', 1];
                    kaf.grow = true;
                end
            end
        end
        
    end
end
