% Test script for kernel adaptive filtering algorithms.
%
% This file tests if the kernel adaptive filtering algorithm is correctly
% implemented.
%
% Instructions:
% 1. Write the name of the algorithm in "algoname". E.g.:
%    >> algoname = 'aldkrls';
% 2. Write a structure containing all the options of the algorithm. E.g.:
%    >> options = struct('nu',.1,'kerneltype','gauss','kernelpar',.2);
% 3. Run this test:
%    >> unit_test;

clearvars -except algoname options
close all
clc

fprintf('Constructing %s object...\n',algoname)
kaf = feval(algoname,options);

fprintf('Testing training phase')
c = 5;
N = 1000;
x = rand(N,2)*c;
y = sin(3*x(:,1)).*cos(x(:,1)+x(:,2));
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end
    kaf = kaf.train(x(i,:),y(i));
end
fprintf('\n')

fprintf('Testing evaluation phase...\n')
[x1,x2] = meshgrid(0:.1:c, 0:.1:c);
yt = kaf.evaluate([x1(:) x2(:)]);
z = reshape(yt,size(x1,1),size(x2,1));
figure;
surf(x1,x2,z);
colormap('spring');
view(20,70)

fprintf('All OK.\n\n')