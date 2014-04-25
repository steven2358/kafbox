% Mixture Kernel Least Mean Square algorithm
%
% Mixture Kernel Least Mean Square (MXKLMS) algorithm, as proposed in R.
% Pokharel, S. Seth, and J.C. Principe, "Mixture kernel least mean square,"
% The 2013 International Joint Conference on Neural Networks (IJCNN),
% pp.1-7, 4-9 Aug. 2013. http://dx.doi.org/10.1109/IJCNN.2013.6706867
%
% Comment: implementation includes a maximum dictionary size M
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

classdef mxklms
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        mu = 2; % learning rate for gradient descent update of v
        eta = .8; % learning rate for updating alpha
        M = 1000; % maximum dictionary size
        kerneltype = 'gauss'; % kernel type
        kernelpars = .5:.5:2; % kernel parameters
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        dict = []; % dictionary
        alpha = []; % expansion coefficients matrix. rows: dict, cols: kern
        psi = []; % kernel weights
        v = []; % gate parameters used to calculate psi
    end
    
    methods
        
        function kaf = mxklms(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if strmatch(fn,fieldnames(kaf),'exact'),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function y_est = evaluate(kaf,x) % evaluate the algorithm
            if size(kaf.dict,1)>0
                N = size(x,1);
                K = multikernel_dict(kaf,x);
                a = kaf.alpha*diag(kaf.psi);
                y_est = reshape(K,N,[])*a(:);
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        function kaf = train(kaf,x,y) % train the algorithm
            if (size(kaf.dict,1)<kaf.M), % avoid infinite growth
                if size(kaf.dict,2)==0 % initialize
                    P = length(kaf.kernelpars);
                    kaf.v = 1/P*ones(P,1);
                    kaf.psi = exp(kaf.v)/sum(exp(kaf.v));
                    kaf.alpha = kaf.eta*kaf.psi'*y;
                    kaf.dict = x;
                else
                    y_est = kaf.evaluate(x);
                    err = y - y_est;
                    
                    % update psi
                    k = multikernel_dict(kaf,x);
                    D = size(kaf.dict,1);
                    yk = sum(reshape(k,D,[]).*kaf.alpha,1)';
                    Jgrad = -2*err/(sum(exp(kaf.v)))^2 * ...
                        (sum(exp(kaf.v))*exp(kaf.v).*yk + ...
                        sum(exp(kaf.v).*yk)*exp(kaf.v));
                    kaf.v = kaf.v - kaf.mu*Jgrad;
                    kaf.psi = exp(kaf.v)/sum(exp(kaf.v));
                    
                    % update alpha
                    kaf.alpha = [kaf.alpha; kaf.eta*kaf.psi'*err]; % grow
                    kaf.dict = [kaf.dict; x]; % grow
                end
            end
        end
        
        function K = multikernel_dict(kaf,X) % multikernel for dictionary
            P = length(kaf.kernelpars); % number of distinct kernels
            N = size(X,1);
            D = size(kaf.dict,1);
            K = zeros(N,D,P);
            switch kaf.kerneltype
                case 'gauss' % RBF kernel
                    norms1 = sum(X.^2,2);
                    norms2 = sum(kaf.dict.^2,2);
                    mat1 = repmat(norms1,1,D);
                    mat2 = repmat(norms2',N,1);
                    
                    d2 = mat1 + mat2 - 2*X*kaf.dict'; % distance matrix
                    Kalpha = -1./(2*kaf.kernelpars.^2);
                    for p=1:P,
                        K(:,:,p) = exp(d2*Kalpha(p));
                    end
                otherwise	% default case
                    for p=1:P,
                        k = kernel(X,kaf.dict,kaf.kerneltype,kaf.kernelpars(p));
                        K(:,:,p) = k;
                    end
            end
        end
    end
end
