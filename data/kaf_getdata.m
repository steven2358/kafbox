function [X,Y] = kaf_getdata(data)

raw = load(sprintf('%s.dat',data.name));

if size(raw,2)==1 % time series
    ydim = 1;
    if isfield(data,'horizon')
        horizon = data.horizon;
    else
        horizon = 1;
    end
elseif size(raw,2)==2 % input output system, input is time series
    ydim = 2;
    horizon = 0;
end

N = min(data.N,size(raw,1)) - horizon;
L = data.L;
X = zeros(N,L);

for i = 1:L,
    X(i:N,i) = raw(1:N-i+1,1); % time embedding
end
Y = raw(1+horizon:N+horizon,ydim); % desired output
