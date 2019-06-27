% L2 algorithm with selective update
%
% M. Ohnishi and M. Yukawa, "Online Learning in L2 Space with Multiple
% Gaussian Kernels", Proceedings of the 25th European Signal Processing
% Conference (EUSIPCO), Kos, Greece, 2017.
%
% Contributed by Motoya Ohnishi
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef L2_s < kernel_adaptive_filter
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        eta = .5; % step size
        mu0 = [0.95,0.95,0.95]; % coherence criterion thresholds
        eps = 1E-2; % regularization
        kerneltype = 'gauss'; % kernel type
        kernelpars = [5,1,0.5]; % kernel parameter
        convex = 0.999; % convex combination
        dictmax = 1000;
        plat = 1e-2; % plat criterion
        largest = 0;
        s = 1;
        rho = 0;
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        dict = []; % combine all dictset
        modict = []; % modulus of the dictionary elements
        alpha = []; % expansion coefficients
        paravector = []; % hyperparameter vector
        typevector = []; % kernel1: 1, kernel2: 2, ...
        Gram = [];
    end
    
    methods
        function kaf = L2_s(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if ismember(fn,fieldnames(kaf)),
                        kaf.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        % evaluate the algorithm
        function y_est = evaluate(kaf,x)
            if size(kaf.dict,1)>0
                k = kaf.multikernel(x,kaf.dict,kaf.kerneltype,0,...
                    zeros(size(x,1),1),kaf.paravector);
                y_est = k*kaf.alpha;
            else
                y_est = zeros(size(x,1),1);
            end
        end
        
        % train the algorithm
        function train(kaf,x,y)
            sizeQ=length(kaf.kernelpars); % number of kernels
            
            if size(kaf.dict,2)==0 % initialize
                kaf.dict = x;
                kx_tmp = kaf.multikernel(x,x,kaf.kerneltype,kaf.largest,...
                    kaf.kernelpars(1),kaf.kernelpars(1));
                kaf.modict = sqrt(kx_tmp);
                kaf.alpha = 0;
                kaf.paravector = kaf.kernelpars(1);
                kaf.typevector = 1;
                kaf.Gram=kx_tmp^(-1);
            else
                for qq=1:sizeQ
                    if(size(kaf.dict,1)==0)
                        C=0;
                        kx_tmp = kaf.multikernel(x,x,kaf.kerneltype,...
                            kaf.largest,kaf.kernelpars(qq),...
                            kaf.kernelpars(qq));
                    else
                        k_tmp = kaf.multikernel(kaf.dict,x,...
                            kaf.kerneltype,kaf.largest,...
                            kaf.paravector,kaf.kernelpars(qq));
                        kx_tmp = kaf.multikernel(x,x,kaf.kerneltype,...
                            kaf.largest,kaf.kernelpars(qq),...
                            kaf.kernelpars(qq));
                        C = k_tmp./(sqrt(kx_tmp)*kaf.modict); % coherence
                    end
                    
                    if (max(C) < kaf.mu0(qq)) % coherence criterion
                        if (qq==1) || ((y-kaf.evaluate(x))^2 > ...
                                kaf.plat*kaf.evaluate(x)^2) % plat
                            if (size(kaf.dict,1)<kaf.dictmax)
                                % when added
                                kaf.dict = [kaf.dict; x];
                                kaf.modict = [kaf.modict; sqrt(kx_tmp)];
                                kaf.alpha = [kaf.alpha; 0];
                                kaf.paravector = [kaf.paravector; ...
                                    kaf.kernelpars(qq)];
                                kaf.typevector = [kaf.typevector; qq];
                                kaf.Gram = kaf.multikernel(kaf.dict,...
                                    kaf.dict,kaf.kerneltype,kaf.largest,...
                                    kaf.paravector,kaf.paravector);
                                break; % only one added
                            else
                                % when exceeding dictmax
                                kaf.dict = [kaf.dict(2:end,:); x];
                                kaf.modict = [kaf.modict(2:end);...
                                    sqrt(kx_tmp)];
                                kaf.alpha = [kaf.alpha(2:end); 0];
                                kaf.paravector = [kaf.paravector(2:end);...
                                    kaf.kernelpars(qq)];
                                kaf.typevector = [kaf.typevector(2:end);...
                                    qq];
                                kaf.Gram = kaf.multikernel(kaf.dict,...
                                    kaf.dict,kaf.kerneltype,kaf.largest,...
                                    kaf.paravector,kaf.paravector);
                                break;
                            end
                        end
                    end
                end
            end
            % selected update
            num_d = min(size(kaf.dict,1),kaf.s);
            dict_id = zeros(num_d,1);
            
            k = kaf.multikernel(x,kaf.dict,kaf.kerneltype,0,0,...
                kaf.paravector);
            
            % selected update
            [~,ii] = sort(k./kaf.modict','descend');
            dict_id(:,1) = ii(1:num_d); %
            
            % selected Gram
            selectGram=kaf.Gram;
            selectGram = selectGram(dict_id(:,1),dict_id(:,1));
            selectGram = kaf.convex*selectGram+(1-kaf.convex)*...
                eye(size(selectGram));
            
            kaf.alpha(dict_id(:,1)) = kaf.alpha(dict_id(:,1)) + ...
                kaf.eta / (kaf.eps + k(dict_id(:,1)) * ...
                (selectGram \ k(dict_id(:,1))')) * ...
                sign(y - k*kaf.alpha) * ...
                max(abs(y - k*kaf.alpha)-kaf.rho,0) * ...
                (selectGram \ k(dict_id(:,1))');
        end
    end
    
    methods (Static = true) % helper functions
        % calculate multikernel
        function K = multikernel(X1,X2,ktype,kpar1,kpardict1,kpardict2)
            N1 = size(X1,1);
            N2 = size(X2,1);
            
            switch ktype
                case 'gauss' % RBF kernel
                    norms1 = sum(X1.^2,2);
                    norms2 = sum(X2.^2,2);
                    
                    sigdict1 = kpardict1.^2;
                    sigdict2 = kpardict2.^2;
                    sigdict = repmat(sigdict1,1,N2) + ...
                        repmat(sigdict2',N1,1)-(kpar1^2*ones(N1,N2));
                    
                    mat1 = repmat(norms1,1,N2);
                    mat2 = repmat(norms2',N1,1);
                    
                    dist2 = mat1 + mat2 - 2*X1*X2';	% full distance matrix
                    K = (2*pi*sigdict).^(-0.5*size(X1,2)).*...
                        exp(-dist2./(2*sigdict));
                    
                otherwise % default case
                    error ('unknown kernel type')
            end
        end
    end
end
