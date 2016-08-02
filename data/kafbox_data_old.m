% KAFBOX_DATA Returns input-output data specified in the options. Either
% load a raw data file or synthetically generate a data set.
%
% The single input argument is a structure that contains the options as
% fields. Possible options:
%   - file: if present, its value is the filename of the data to load.
%   - generate: if present, its value is the name of the function used to
%   generate the data.
%   - horizon: prediction horizon
%   - N: number of data to load. Limited by the number of data available
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

function [X,Y,X_test,Y_test] = kafbox_data_old(options)

if isfield(options,'file')
    
    data = load(options.file);
    
    switch class(data)
        case 'double' % raw numeric data
            
            % time embedding
            embedding = 0;
            if isfield(options,'embedding')
                embedding = options.embedding;
            end
            
            if size(data,2)==1 % time series
                
                % prediction horizon
                horizon = 1;
                if isfield(options,'horizon')
                    horizon = options.horizon;
                end
                
                % construct signal
                X = time_embedding(data,embedding);
                Y = data(1+horizon:end); % desired output
                
            elseif size(data,2)==2 % input output system; input is time series
                
                % construct signal
                X = time_embedding(data(:,1),embedding);
                Y = data(1:end,2); % desired output
                
            end
            
            % apply offset
            if isfield(options,'offset')
                X = X(1+options.offset:end,:);
                Y = Y(1+options.offset:end);
            end
            
            % crop
            N = length(Y);
            if isfield(options,'N')
                N = min(options.N,N);
            end
            X = X(1:N,:);
            Y = Y(1:N);
            
            X_test = [];
            Y_test = [];
            
        case 'struct'
            
            X = data.X_train;
            Y = data.Y_train;
            X_test = data.X_test;
            Y_test=  data.Y_test;
            
            if isfield(options,'permutation')
                if options.permutation ~= 0,
                    
                    randseed = options.permutation;
                    randn('state',randseed); %#ok<RAND>
                    rand('state',randseed); %#ok<RAND>
                    
                    N = size(X,1);
                    rp = randperm(N);
                    X = X(rp,:);
                    Y = Y(rp);
                    
                end
            end
            
            % crop
            N = length(Y);
            if isfield(options,'N')
                N = min(options.N,N);
            end
            X = X(1:N,:);
            Y = Y(1:N);
            
        otherwise
            error('unknown data type');
    end
elseif isfield(options,'generate')
    
    if isfield(options,'N_test')
        eval(sprintf('[X,Y,y_ref,X_test,Y_test] = generate_%s(options);',...
            options.generate));
    else
        eval(sprintf('[X,Y] = generate_%s(options);',options.generate));
        X_test = [];
        Y_test = [];
    end
    
end


function X_embedded = time_embedding(X,L)
N = size(X,1);
X_embedded = zeros(N,L);
for i = 1:L,
    X_embedded(i:N,i) = X(1:N-i+1,1); % time embedding
end
