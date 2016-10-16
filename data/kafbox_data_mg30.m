% KAFBOX_DATA_MG30 Data loader for Mackey-Glass 30 data set

function [X,Y,X_test,Y_test] = kafbox_data_mg30(data_opts)

data = load('mg30.dat');

% prediction horizon
horizon = 1;
if isfield(data_opts,'horizon')
    horizon = data_opts.horizon;
end

% construct signal
X = data(1:end-horizon); % input
Y = data(1+horizon:end); % desired output

X_test = [];
Y_test = [];
