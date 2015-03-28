% Kernel Recursive Least-Squares Tracker algorithm
%
% S. Van Vaerenbergh, M. Lazaro-Gredilla, and I. Santamaria, "Kernel
% Recursive Least-Squares Tracker for Time-Varying Regression," IEEE
% Transactions on Neural Networks and Learning Systems, vol. 23, no. 8, pp.
% 1313-1326, Aug. 2012, http://dx.doi.org/10.1109/TNNLS.2012.2200500
%
% Comment: using back-to-the-prior forgetting
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef krlst
    
    properties (GetAccess = 'public', SetAccess = 'private')
        lambda = .999; % forgetting factor
        sn2 = 1E-2; % noise to signal ratio
        M = 100; % dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
        jitter = 1E-6; % jitter noise to avoid roundoff error
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        dict = []; % dictionary
        Q = []; % inverse kernel matrix
        mu = []; % posterior mean
        Sigma = []; % posterior covariance
        nums02ML = 0;
        dens02ML = 0;
        s02 = 0; % signal power, adaptively estimated
        prune = false; % flag
        reduced = false; % flag
    end
    
    methods
        
        function kaf = krlst(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if strmatch(fn,fieldnames(kaf),'exact'),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function [mean_test, var_test] = evaluate(kaf,x) % evaluate
            if size(kaf.dict,1)>0
                
                k = kernel(kaf.dict,x,kaf.kerneltype,kaf.kernelpar);
                q = kaf.Q*k;
                mean_test = q'*kaf.mu; % predictive mean
                
                if nargout>1
                    ktt = kernel(x,x,[kaf.kerneltype '-diag'],kaf.kernelpar);
                    sf2 = ktt + kaf.jitter + sum(k.*((kaf.Q*kaf.Sigma*kaf.Q-kaf.Q)*k),1)';
                    sf2(sf2<0) = 0;
                    var_test = kaf.s02*(kaf.sn2 + sf2); % predictive variance
                end
            else
                mean_test = zeros(size(x,1),1); % prior
                % var_test % prior
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            m = size(kaf.Sigma,1);
            
            if m<1 % initialize
                k = kernel(x,x,kaf.kerneltype,kaf.kernelpar);
                k = k + kaf.jitter;
                kaf.Q = 1/k;
                kaf.mu = y*k/(k+kaf.sn2);
                kaf.Sigma = k - k^2/(k+kaf.sn2);
                kaf.dict = x; % dictionary bases
                kaf.nums02ML = y^2/(k+kaf.sn2);
                kaf.dens02ML = 1;
                kaf.s02 = kaf.nums02ML / kaf.dens02ML;
            else
                % forget
                K = kernel(kaf.dict,kaf.dict,kaf.kerneltype,kaf.kernelpar);
                kaf.Sigma = kaf.lambda*kaf.Sigma + (1-kaf.lambda)*K; % forget
                kaf.mu = sqrt(kaf.lambda)*kaf.mu; % forget
                
                % predict
                k = kernel([kaf.dict; x],x,kaf.kerneltype,kaf.kernelpar);
                kt = k(1:end-1);
                ktt = k(end) + kaf.jitter;
                q = kaf.Q*kt;
                y_mean = q'*kaf.mu; % predictive mean
                gamma2 = ktt - kt'*q; gamma2(gamma2<0)=0; % projection uncertainty
                h = kaf.Sigma*q;
                sf2 = gamma2 + q'*h; sf2(sf2<0)=0; % noiseless prediction variance
                sy2 = kaf.sn2 + sf2;
                % y_var = s02*sy2; % predictive variance
                
                % include new sample and add a basis
                Qold = kaf.Q;
                p = [q; -1];
                kaf.Q = [kaf.Q zeros(m,1);zeros(1,m) 0] + 1/gamma2*(p*p');
                
                p = [h; sf2];
                kaf.mu = [kaf.mu; y_mean] + (y - y_mean)/sy2*p; % posterior mean
                kaf.Sigma = [kaf.Sigma h; h' sf2] - 1/sy2*(p*p'); % posterior covariance
                m = m + 1;
                kaf.dict = [kaf.dict; x];
                
                % estimate s02 via ML
                kaf.nums02ML = kaf.nums02ML + kaf.lambda*(y - y_mean)^2/sy2;
                kaf.dens02ML = kaf.dens02ML + kaf.lambda;
                kaf.s02 = kaf.nums02ML/kaf.dens02ML;
                
                kaf.prune = false;
                % delete a basis if necessary
                if (m>kaf.M  || gamma2<kaf.jitter)
                    if gamma2<kaf.jitter, % to avoid roundoff error
                        if gamma2<kaf.jitter/10
                            warning('Numerical roundoff error too high, you should increase jitter noise') %#ok<WNTAG>
                        end
                        criterion = [ones(1,m-1) 0];
                    else % MSE pruning criterion
                        errors = (kaf.Q*kaf.mu)./diag(kaf.Q);
                        criterion = abs(errors);
                    end
                    [~, r] = min(criterion); % remove element r, which incurs in the minimum error
                    smaller = 1:m; smaller(r) = [];
                    
                    if r == m, % if we must remove the element we just added, perform reduced update instead
                        kaf.Q = Qold;
                        kaf.reduced = true;
                    else
                        Qs = kaf.Q(smaller, r);
                        qs = kaf.Q(r,r); kaf.Q = kaf.Q(smaller, smaller);
                        kaf.Q = kaf.Q - (Qs*Qs')/qs;
                        kaf.reduced = false;
                    end
                    kaf.mu = kaf.mu(smaller);
                    kaf.Sigma = kaf.Sigma(smaller, smaller);
                    kaf.dict = kaf.dict(smaller,:);
                    kaf.prune = true;
                end
            end
        end
        
    end
end
