% Wobbly field regression demo. Also used in lib/test/unit_test.m
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox

close all
clear

%% PARAMETERS

algorithm = 'krlst'; % algorithm class (choose from lib/ folder)
opts = struct(); % algorithm options go here (kernel type, parameters, etc)

%% PROGRAM

kaf = feval(algorithm, opts);

% generate some data
c = 5;
N = 1000;
x = rand(N,2)*c;
y = sin(3*x(:,1)).*cos(x(:,1)+x(:,2));

fprintf('Training')
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end
    kaf = kaf.train(x(i,:),y(i));
    % y_test = kaf.evaluate(x(i+1,:));
end
fprintf('\n')

%% OUTPUT

[x1,x2] = meshgrid(0:.2:c, 0:.2:c);
yt = kaf.evaluate([x1(:) x2(:)]);

z = reshape(yt,size(x1,1),size(x2,1));
figure;
surf(x1,x2,z);
