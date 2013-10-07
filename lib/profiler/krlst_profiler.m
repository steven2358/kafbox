% Profiler extension for Kernel Recursive Least-Squares Tracker algorithm
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef krlst_profiler < krlst
    
    properties (GetAccess = 'public', SetAccess = 'private')
        elapsed = 0; % elapsed time
    end
    
    methods
        
        function kaf = krlst_profiler(parameters) % constructor
            if nargin<1, parameters = struct(); end
            kaf = kaf@krlst(parameters);
        end
        
        function flops = lastflops(kaf) % flops for last iteration
            m = size(kaf.dict,1);
            if kaf.prune
                m1 = m;
            else
                m1 = m - 1;
            end
            if kaf.reduced
                m3 = 0;
            else
                m3 = m;
            end
            m2 = m1 + 1;
            
            floptions = struct(...
                'sum', m1^2 + m1 + 1 + m1^2 - m1 + m1 - 1 + m1 + m1 - 1 + m1 + 1 + m2^2 + m2 + 1 + m2^2 + 2 + 1 + m2^2 - m2 + m3^2, ...
                'mult', 2*m1^2 + m1^2 + m1 + m1 + m1 + m1 + m2^2 + 1 + 1 + m2^2 + 3 + m2^2 + m3^2, ...
                'div', 1 + 1 + 2 + m2 + m3, ...
                sprintf('%s_kernel',kaf.kerneltype), [m1^2 + 1,1,size(kaf.dict,2)]);
            
            flops = kflops(floptions);
        end
        
        %% flops breakdown
        
        % K = kernel(kaf.dict,kaf.dict,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1^2
        
        % kaf.Sigma = kaf.lambda*kaf.Sigma + (1-kaf.lambda)*K + kaf.jitter*eye(m); % forget
        % sum: m1^2 + m1 + 1
        % mult: 2*m1^2
        
        % kaf.mu = sqrt(kaf.lambda)*kaf.mu; % square can be pre-calculated
        % mult: m1
        
        % k = kernel([kaf.dict; x],x,kaf.kerneltype,kaf.kernelpar);
        % kernel: m1
        
        % q = kaf.Q*kt;
        % sum: m1^2 - m1
        % mult: m1^2
        
        % y_mean = q'*kaf.mu; % predictive mean
        % sum: m1 - 1
        % mult: m1
        
        % gamma2 = ktt - kt'*q; gamma2(gamma2<0)=0; % projection uncertainty
        % sum: m1
        % mult: m1
        
        % h = kaf.Sigma*q;
        % sum: m1 - 1
        % mult: m1
        
        % sf2 = gamma2 + q'*h; sf2(sf2<0)=0; % noiseless prediction variance
        % sum: m1
        % mult: m1
        
        % sy2 = kaf.sn2 + sf2;
        % sum: 1
        
        % kaf.Q = [kaf.Q zeros(m,1);zeros(1,m) 0] + 1/gamma2*(p*p'); with p = [q; -1];
        % sum: m2^2
        % mult: m2^2 + 1
        % div: 1
        
        % kaf.mu = [kaf.mu; y_mean] + (y - y_mean)/sy2*p; % with p = [h; sf2];
        % sum: m2 + 1
        % mult: m2 + 1
        % div: 1
        
        % kaf.Sigma = [kaf.Sigma h; h' sf2] - 1/sy2*(p*p'); % posterior covariance
        % sum: m2^2
        % mult: m2^2
        
        % kaf.nums02ML = kaf.nums02ML + kaf.lambda*(y - y_mean)^2/sy2;
        % sum: 2
        % mult: 3
        
        % kaf.dens02ML = kaf.dens02ML + kaf.lambda;
        % sum: 1
        
        % kaf.s02 = kaf.nums02ML/kaf.dens02ML;
        % div: 1
        
        % errors = (kaf.Q*kaf.mu)./diag(kaf.Q); % MSE pruning criterion
        % sum: m2^2 - m2
        % mult: m2^2
        % div: m2
        
        % kaf.Q = kaf.Q - (Qs*Qs')/qs; % if removed element is not the last
        % sum: m3^2
        % mult: m3^2
        % div: m3
        
        %%
        
        function kaf = train_profiled(kaf,x,y)
            % kaf.prev_dict_size = size(kaf.dict,1);
            t1 = tic;
            kaf = kaf.train(x,y);
            t2 = toc(t1);
            kaf.elapsed = kaf.elapsed + t2;
        end
        
        function bytes = lastbytes(kaf) % bytes used in last iteration
            m = size(kaf.dict,1);
            bytes = 8*(m^2 + m^2 + m + 2 + m*size(kaf.dict,2)); % 8 bytes for double precision
            % Q, Sigma, mu, nums02ML, dens02ML, dict
        end
        
    end
end
