% Showcasing Extended kernel RLS
%
% "Problem" with this implementation is that is does not update anymore
% once number of presented data points exceed M.

close all
clear

%% PARAMETERS

algorithm = 'exkrls'; % algorithm class (choose from lib/ folder)
opts = struct(
  "alphaf", .999, % state forgetting factor, "alpha" in publication
  "beta", .995, % data forgetting factor
  "lambda", 1E-2, % regularization parameter
  "q", 1E-3, % trade-off between modeling variation and measurement disturbance
  "M", 100,  # dictionary size  # not necessarily the lookback
  # "kerneltype", 'gauss'
  "kerneltype", 'linear'  # can only fit a plane

); % algorithm options go here (kernel type, parameters, etc)

%% PROGRAM

kaf = feval(algorithm, opts); %#ok<FVAL>

% generate some data
c = 5;
N = 200;
x = rand(N,2)*c;
y = zeros(N,1);
% y = sin(3*x(:,1)).*cos(x(:,1)+x(:,2));  # non-linear surface
y(1:N/2) = 3*x(1:N/2,1)+0.5*x(1:N/2,2);  # linear surface
y(N/2:N) = -1*x(N/2:N,1)-1.5*x(N/2:N,2);  # linear surface changes
y = y + randn(N,1)*1.0;  # add some noise

fprintf('First half of Training')
for i=1:(N/2-1)
    if ~mod(i,floor(N/10)), fprintf('.'); end
    kaf.train(x(i,:),y(i));
    % y_test = kaf.evaluate(x(i+1,:));
end
fprintf('\n')

%% OUTPUT

[x1,x2] = meshgrid(0:.2:c, 0:.2:c);
yt = kaf.evaluate([x1(:) x2(:)]);

z = reshape(yt,size(x1,1),size(x2,1));
figure;
surf(x1,x2,z);
hold;
plot3(x(:,1), x(:,2), y, "+");

fprintf('Training shortly after Non-Stationarity')
for i=N/2:N/2+9
    if ~mod(i,floor(N/10)), fprintf('.'); end
    kaf.train(x(i,:),y(i));
    % y_test = kaf.evaluate(x(i+1,:));
end
fprintf('\n')

%% OUTPUT

[x1,x2] = meshgrid(0:.2:c, 0:.2:c);
yt = kaf.evaluate([x1(:) x2(:)]);

z = reshape(yt,size(x1,1),size(x2,1));
figure;
surf(x1,x2,z);
hold;
plot3(x(:,1), x(:,2), y, "+");

fprintf('Training Rest')
for i=N/2+10:N
    if ~mod(i,floor(N/10)), fprintf('.'); end
    kaf.train(x(i,:),y(i));
    % y_test = kaf.evaluate(x(i+1,:));
end
fprintf('\n')

%% OUTPUT

[x1,x2] = meshgrid(0:.2:c, 0:.2:c);
yt = kaf.evaluate([x1(:) x2(:)]);

z = reshape(yt,size(x1,1),size(x2,1));
figure;
surf(x1,x2,z);
hold;
plot3(x(:,1), x(:,2), y, "+");
