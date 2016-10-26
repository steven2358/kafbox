% KAFBOX_DATA Data handler. Returns input-output data specified in the
% options.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

function [X_train,Y_train,X_test,Y_test] = kafbox_data(data_options)

if isfield(data_options,'class')
    data_name = lower(data_options.class);
elseif isfield(data_options,'name')
    data_name = lower(data_options.name);
else
    error('Unknown data.');
end

dataset_handle = str2func(['kafbox_data_' data_name]);

[X,Y,X_test,Y_test] = feval(dataset_handle, data_options);

% time embedding
L = 0;
if isfield(data_options,'embedding')
    L = data_options.embedding;
    X = time_embedding(X,L);
end

% apply offset
if isfield(data_options,'offset')
    X = X(1+data_options.offset:end,:);
    Y = Y(1+data_options.offset:end);
end

% crop training data
N = length(Y);
if isfield(data_options,'N')
    N = min(data_options.N,N);
end
X_train = X(1:N,:);
Y_train = Y(1:N);

if isfield(data_options,'n_train') && isfield(data_options,'n_test')
    n_train = data_options.n_train;
    n_test = data_options.n_test;
    X_train = X(L:n_train+L-1,:);
    Y_train = Y(L:n_train+L-1);
    X_test = X(n_train+L:n_train+L+n_test-1,:);
    Y_test = Y(n_train+L:n_train+L+n_test-1);
end


function X_embedded = time_embedding(X,L)
N = size(X,1);
X_embedded = zeros(N,L); % L-1 first vectors are zero padded
for i = 1:L,
    X_embedded(i:N,i) = X(1:N-i+1,1); % time embedding
end
