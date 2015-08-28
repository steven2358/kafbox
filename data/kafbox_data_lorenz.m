% KAFBOX_DATA_LORENZ Data loader for Lorenz data set

function [X,Y,X_test,Y_test] = kafbox_data_lorenz(data_opts)

data = load('lorenz.dat');

% prediction horizon
horizon = 1;
if isfield(data_opts,'horizon')
    horizon = data_opts.horizon;
end

% construct signal
X = data(horizon:end-1); % input
Y = data(1+horizon:end); % desired output

X_test = [];
Y_test = [];
