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
end

