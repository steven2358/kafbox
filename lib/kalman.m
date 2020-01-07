% Discrete-time Kalman filter.
%
% R. E. Kalman, "A New Approach to Linear Filtering and Prediction
% Problems", ASME Journal of Basic Engineering, 1960, 82(1), pp. 35-45.
% http://dx.doi.org/10.1115/1.3662552
%
% Notation is taken from: R. Faragher, "Understanding the Basis of the
% Kalman Filter Via a Simple and Intuitive Derivation [Lecture Notes]," in
% IEEE Signal Processing Magazine, vol. 29, no. 5, pp. 128-132, Sept. 2012.
% http://dx.doi.org10.1109/MSP.2012.2203621
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

classdef kalman < linear_filter
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        F; % state transition matrix
        B; % control input matrix
        Q; % state noise covariance matrix
        % H; % state-to-measurements transition matrix
        R; % measurement noise covariance matrix
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        x_mean; % state mean
        P; % state covariance
    end
    
    methods
        function obj = kalman(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)',
                    if ismember(fn,fieldnames(obj)),
                        obj.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function [mean_test, var_test] = evaluate(obj,H) % evaluate
            if ~isempty(obj.x_mean)
                mean_test = H*obj.x_mean;
                var_test = H*obj.P*H'; % predictive variance
            else
                mean_test = zeros(size(H,1),1);
                var_test = nan(size(H,1),1); % signal scale is unknown
            end
        end
        
        function train(obj,H,z,u) % train the algorithm
            if isempty(obj.x_mean) % initialize
                [n,m] = size(H); % measurement and state dimensions
                obj.x_mean = zeros(m,1); % filter coefficients
                obj.P = eye(m);
                
                if isempty(obj.F)
                    % default state transition: random walk
                    obj.F = eye(m);
                    obj.Q = 0.1*eye(m);
                    % default control input
                    obj.B = eye(m);
                    % default measurement noise covariance
                    obj.R = 1E-2*eye(n);
                end
            end
            
            if nargin<4
                u = zeros(size(H,2),1); % no control input
            end
            
            % state prediction equations
            obj.x_mean = obj.F*obj.x_mean + obj.B*u;
            obj.P = obj.F*obj.P*obj.F + obj.Q;
            
            % calculate the kalman gain
            K = obj.P*H'/(H*obj.P*H' + obj.R);
            
            % measurement update equations
            obj.x_mean = obj.x_mean + K*(z - H*obj.x_mean); % mean
            obj.P = obj.P - K*H*obj.P; % covariance
        end
    end
end
