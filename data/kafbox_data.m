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
% http://sourceforge.net/projects/kafbox/

function [X,Y] = kafbox_data(options)

if isfield(options,'file')
    
    raw = load(options.file);
    
    if size(raw,2)==1 % time series
        ydim = 1;
        if isfield(options,'horizon')
            horizon = options.horizon;
        else
            horizon = 1;
        end
    elseif size(raw,2)==2 % input output system, input is time series
        ydim = 2;
        horizon = 0;
    end
    
    if isfield(options,'N')
        N = min(options.N,size(raw,1)) - horizon;
    else
        N = size(raw,1) - horizon;
    end
    
    L = options.embedding;
    X = zeros(N,L);
    for i = 1:L,
        X(i:N,i) = raw(1:N-i+1,1); % time embedding
    end
    Y = raw(1+horizon:N+horizon,ydim); % desired output
   
elseif isfield(options,'generate')
    
    eval(sprintf('[X,Y] = generate_%s(options);',options.generate));
    
end


