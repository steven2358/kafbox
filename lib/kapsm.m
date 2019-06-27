% Kernel Affine Projection Subgradient Method
%
% K. Slavakis, S. Theodoridis, and I. Yamada, "Online kernel-based
% classification using adaptive projection algorithms," IEEE Transactions
% on Signal Processing, Vol. 56, No. 7, pp. 2781-2796, 2008.
% http://dx.doi.org/10.1109/TSP.2008.917376
%
% Remark: implemented with L2-ball forgetting. Code contributed by
% Pantelis Bouboulis.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef kapsm < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        M = 200; % dictionary size
        Delta = 10; % L2-ball forgetting radius
        mu = 1.8, % learning rate
        Q = 10; % number of subgradients
        loss = 'l2'; % loss function type
        loss_param = 2; % Huber loss parameter
        kerneltype = 'gauss'; % kernel type
        kernelpar = .5; % kernel parameter
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        dict = []; % dictionary
        alpha = []; % expansion coefficients
        norm_f = []; % function norm for L2-ball forgetting
        xmem = []; % input memory
        ymem = []; % output memory
    end
    
    methods
        function kaf = kapsm(parameters) % constructor
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
            if size(kaf.xmem,1)<kaf.Q
                % grow memory
                kaf.xmem = [kaf.xmem; x];
                kaf.ymem = [kaf.ymem; y];
            else
                % slide memory
                kaf.xmem = [kaf.xmem(2:kaf.Q,:); x];
                kaf.ymem = [kaf.ymem(2:kaf.Q); y];
            end
            q = size(kaf.xmem,1);
            
            d_hat = kaf.evaluate(kaf.xmem);
            e = kaf.ymem - d_hat;
            
            L_n = kaf.loss_fun(e,kaf.loss,kaf.loss_param);
            subgrad = kaf.compute_subgrad_coef(e,kaf.loss,kaf.loss_param);
            norm_subgrad = subgrad.^2 .* ...
                kernel(kaf.xmem(end-q+1:end,:),...
                kaf.xmem(end-q+1:end,:),...
                [kaf.kerneltype '-diag'],kaf.kernelpar);
            
            omega = 1/q*ones(q,1); % convex weights
            
            beta = omega.*L_n.*subgrad./(norm_subgrad+eps);
            
            % % extrapolation parameter
            % nominator = (omega.*L_n)'*(L_n./(norm_subgrad+eps));
            % K = kernel(kaf.xmem,kaf.xmem,kaf.kerneltype,kaf.kernelpar);
            % denominator = beta'*K*beta;
            % Mn = nominator/(denominator+eps);
            % mu = .5*Mn;
            
            kaf.alpha = [kaf.alpha; 0];
            kaf.alpha(end-q+1:end) = kaf.alpha(end-q+1:end) - kaf.mu*beta;
            kaf.dict = [kaf.dict; x];
            
            % L2-ball forgetting: update function norm
            K = kernel(kaf.xmem,kaf.xmem,kaf.kerneltype,kaf.kernelpar);
            R_sum = beta(end-q+1:end)'*K*beta(end-q+1:end);
            sum2 = beta(1:end-q)'*d_hat(1:end-q);
            kaf.norm_f = sqrt(kaf.norm_f*kaf.norm_f + 2*sum2 + R_sum);
            if kaf.norm_f > kaf.Delta
                kaf.alpha = kaf.Delta/kaf.norm_f * kaf.alpha;
                kaf.norm_f = kaf.Delta;
            end
            
            % sliding-window pruning
            if length(kaf.alpha)>kaf.M,
                k = kernel(kaf.dict,kaf.dict(1,:),...
                    kaf.kerneltype,kaf.kernelpar);
                kaf.norm_f = sqrt(kaf.norm_f.^2 + kaf.alpha(1)^2*k(1) - ...
                    2*kaf.alpha(1)*k'*kaf.alpha);
                kaf.dict(1,:) = [];
                kaf.alpha(1) = [];
            end
        end
        
    end
    
    methods (Static = true)
        
        function loss = loss_fun(ksi,loss_type,loss_param)
            switch loss_type
                case 'l2'
                    loss = ksi.^2;
                case 'l1'
                    loss = abs(ksi);
                case 'huber'
                    sigma = loss_param;
                    loss = zeros(length(ksi),1);
                    for i=1:length(ksi),
                        if abs(ksi(i)) <= sigma
                            loss(i) = ksi(i)^2/(2*sigma);
                        else
                            loss(i) = abs(ksi(i)) - sigma/2;
                        end
                    end
            end
        end
        
        function subgrad_coef = ...
                compute_subgrad_coef(e,loss_type,loss_param)
            subgrad_coef = zeros(length(e),1);
            for i=1:length(e)
                switch loss_type
                    case 'l2'
                        subgrad_coef(i) = -2*e(i);
                    case 'l1'
                        subgrad_coef(i) = -sign(e(i));
                    case 'huber'
                        sigma = loss_param;
                        if abs(e(i)) <= sigma
                            subgrad_coef(i) = -e(i)/sigma;
                        else
                            subgrad_coef(i) = -sign(e(i));
                        end
                end
            end
        end
    end
end
