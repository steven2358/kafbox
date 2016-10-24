% Parallel HYperslab Projection along Affine SubSpace (PHYPASS) algorithm
%
% M. Takizawa and M. Yukawa, "An Efficient Data-Reusing Kernel Adaptive
% Filtering Algorithm Based on Parallel Hyperslab Projection Along Affine
% Subspace," in Proc. ICASSP, May. 2013, pp.3557-3561.
%
% Remark: version01, August 2013
% Contributor for this code: Masa-aki Takizawa
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef phypass < handle
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mu = 0.5; % step size
        s = 1; % number of update coefficients (dictionary size)
        sigma = 0.95; % threshold for dictionary
        p = 8; % number of hyper slabs
        kerneltype = 'gauss'; % kernel type
        kernelpar = 1; % kernel parameter
        omega = 1/8; % weight
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
        mem = []; % memory
        dict = []; % dictionary
        alpha = []; % expansion coefficients
        d = []; % output window
    end
    
    methods
        
        function kaf = phypass(parameters) % constructor
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
            if (length(kaf.d) < kaf.p) % grow the memory
                kaf.mem = [kaf.mem; x];
                kaf.d = [kaf.d; y];
            else % sliding window
                kaf.mem = [kaf.mem(2:end,:); x];
                kaf.d = [kaf.d(2:end); y];
            end
            if size(kaf.dict,1) == 0 % initialize
                kaf.dict = x;
                kaf.alpha = 0;
            else
                k = kernel(x,kaf.dict,kaf.kerneltype,kaf.kernelpar);
                if max(k) < kaf.sigma % coherence criterion
                    kaf.dict = [kaf.dict ; x]; % add the current input x
                    kaf.alpha = [kaf.alpha ; 0];
                end
            end
            num = length(kaf.d); % length of memory
            
            num_d = min(size(kaf.dict,1),kaf.s); % number of dictionary elements to update
            
            dict_id = zeros(num_d,num); % memory of dictionary index
            
            Kmemdict = kernel(kaf.mem,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            Kdict = kernel(kaf.dict,kaf.dict,kaf.kerneltype,kaf.kernelpar);
            
            for k=1:num,
                d_check = Kmemdict(k,:); % kernel between k'th memory element and full dictionary
                [mm,ii] = sort(d_check,'descend'); %#ok<ASGLU>
                dict_id(:,k) = ii(1:num_d); % memory indexes
            end
            
            k_n = zeros(num_d,num_d,num);
            y_n = zeros(num_d,num);
            alpha_new = zeros(num_d,num);
            for k=1:num
                k_n(:,:,k) = Kdict(dict_id(:,k),dict_id(:,k));
                y_n(:,k) = Kmemdict(k,dict_id(:,k))';
                alpha_new(:,k) = k_n(:,:,k)\y_n(:,k); % compute projections onto dictionary subspaces
            end
            ln_nume = 0;
            numerator = zeros(num,1);
            denominator = zeros(num,1);
            beta = zeros(num,1);
            
            for k=1:num,
                numerator(k) = Kmemdict(k,:)*kaf.alpha;
                denominator(k) = Kmemdict(k,dict_id(:,k))*alpha_new(:,k);
                beta(k) = (kaf.d(k) - numerator(k))/(denominator(k)); %  progress of the projection onto hyperplanes
                ln_nume = ln_nume + kaf.omega * beta(k)^2 *...
                    alpha_new(:,k)' * k_n(:,:,k) * alpha_new(:,k); % numerator of the extrapolation coefficcient
            end
            
            G = zeros(num,num);
            for k = 1:num,
                Grow =  alpha_new(:,k)'*Kdict(dict_id(:,k),dict_id);
                sti = 1;
                ndi = num_d;
                for l = 1:num,
                    G(k,l) = Grow(1,sti:ndi)*alpha_new(:,l);
                    sti = ndi + 1;
                    ndi = ndi + num_d;
                end
            end
            
            ln_deno = kaf.omega^2 * beta' * G * beta; % denominator of the extrapolation coefficient
            ln = ln_nume/(ln_deno); % extrapolation coefficient
            alpha_new2 = zeros(length(kaf.alpha),num);
            for k=1:num
                alpha_new2(dict_id(:,k),k) = alpha_new(:,k);
            end
            for k=1:num
                kaf.alpha = kaf.alpha + kaf.mu * ln *...
                    kaf.omega * beta(k) * alpha_new2(:,k); % update expansion coefficients
            end
            
        end
    end
end
