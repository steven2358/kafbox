% Multilayer Perceptron.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/
%
% Author: Steven Van Vaerenbergh, 2006.

classdef mlp < base_estimator
    
    properties (GetAccess = 'public', SetAccess = 'private') % parameters
        eta = 0.01; %learning rate
        n_input; % number of inputs
        n_hidden = 10; % neurons in hidden layer
        output_scale = 1; % scaling factor for output
    end
    
    properties (GetAccess = 'public', SetAccess = 'private') % variables
        b1; % bias for hidden layer
        W1; % weights for hidden layer
        b2; % bias for output layer
        W2; % weights for output layer
    end
    
    methods
        function obj = mlp(parameters) % constructor
            if (nargin > 0) % copy valid parameters
                for fn = fieldnames(parameters)'
                    if ismember(fn,fieldnames(obj))
                        obj.(fn{1}) = parameters.(fn{1});
                    end
                end
            end
        end
        
        function [y_est,y1,y1p] = evaluate(obj,xT) % evaluate the algorithm
            if ~isempty(obj.W1)
                z = obj.W1*xT' + repmat(obj.b1,1,size(xT,1));
                [y1,y1p] = obj.nl(z);
                y_est = (y1'*obj.W2 + obj.b2')*obj.output_scale;
            else
                % not initialized
                y_est = zeros(size(xT,1),1);
                yl = zeros(size(xT,1),1);
                ylp = zeros(size(xT,1),1);
            end
        end
        
        function initialize(obj,n_input)
            obj.n_input = n_input;
            % Xavier initialization
            ni = obj.n_input;
            no = obj.n_hidden;
            obj.W1 = 1/sqrt(ni+no)*randn(no,ni);
            obj.b1 = zeros(no,1);
            obj.W2 = 1/sqrt(no+1)*randn(no,1);
            obj.b2 = 0;
        end
        
        function train(obj,xT,y) % train the algorithm
            if isempty(obj.W1)
                obj.initialize(length(xT));
            end
            x = xT';
            [y_est,y1,y1p] = obj.evaluate(xT);
            
            e = (y - y_est)/obj.output_scale; % linear output
            [dVdW1,dVdW2,dVdb1,dVdb2] = mlp.dvdw_mse(x,y1p,y1,e,obj.W2);
            
            % evaluate gradient at current weight vector
            obj.W1 = obj.W1 - obj.eta*dVdW1;
            obj.W2 = obj.W2 - obj.eta*dVdW2;
            obj.b1 = obj.b1 - obj.eta*dVdb1;
            obj.b2 = obj.b2 - obj.eta*dVdb2;
        end
    end
    
    methods (Static = true)
        
        function [y,yp] = nl(x)
            % y = f(x), yp = f'(x)
            
            y = tanh(x);yp = ones(size(y))-y.^2; % d(tanh(x))/dx=1-tanh(x)^2
        end
        
        function [dy2dW1,dy2db1,dy2dW2,dy2db2] = dydw(x,y1p,y1,W2)
            
            N = length(x(1,:)); % number of data points in training set
            n0 = length(x(:,1));
            
            dy2db2 = ones(1,N);
            dy2dW2 = zeros(size(y1));
            dy2db1 = zeros(size(y1));
            dy2dW1 = zeros(size(y1,1),n0,size(y1,2));
            for n = 1:N
                dy2dW2(:,n) = y1(:,n);
                dy2db1(:,n) = W2.*y1p(:,n);
                for j = 1:n0
                    dy2dW1(:,j,n) = dy2db1(:,n)*x(j,n);
                end
            end
        end
        
        function [dVdW1,dVdW2,dVdb1,dVdb2] = dvdw_mse(x,y1p,y1,e,W2)
            
            % obtain the derivative of NN output wrt weights
            [dy2dW1,dy2db1,dy2dW2] = mlp.dydw(x,y1p,y1,W2);
            
            N = length(x(1,:)); % number of data points in training set
            n0 = length(x(:,1)); % length of input vector
            n = length(W2); % number of hidden neurons
            
            e = e(:)'; % e must be a row vector
            
            dVdW2 = zeros(n,1);dVdW1 = zeros(n,n0);dVdb1 = zeros(n,1);dVdb2=0;
            for n = 1:N
                dVdW2 = dVdW2 -2* e(n)*dy2dW2(:,n);
                dVdW1 = dVdW1 -2* e(n)*dy2dW1(:,:,n);
                dVdb1 = dVdb1 -2* e(n)*dy2db1(:,n);
                dVdb2 = dVdb2 -2* e(n);
            end
        end
        
    end
end
