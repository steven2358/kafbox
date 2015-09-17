% Partly reproduces figure 2 from "Online Prediction of Time Series Data
% With Kernels". (5 algorithms, 20 MC simulations.)
%
% Learning curves for KAP, KNLMS, NORMA and KRLS on a nonlinear system.
% Execution time: 3.5 minutes (Intel Pentium Core2 Duo).
%
% C. Richard, J.C.M. Bermudez, and P. Honeine, "Online Prediction of Time
% Series Data With Kernels," IEEE Transactions on Signal Processing,
% vol. 57, no. 3, pp. 1058-1067, March 2009,
% http://dx.doi.org/10.1109/TSP.2008.2009895
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab
% http://sourceforge.net/projects/kafbox/

clear
close all

%% PARAMETERS

N = 3000; % number of training data points
ktype = 'laplace';
kpar = 0.35;

setups{1} = norma(struct('lambda',0.09,'tau',35,'eta',0.09,'tcoeff',-1/2,'kerneltype',ktype,'kernelpar',kpar));
setups{2} = knlms(struct('eta',.01,'eps',9E-4,'mu0',0.3,'kerneltype',ktype,'kernelpar',kpar));
setups{3} = kap(struct('p',2,'eta',.009,'eps',7e-2,'mu0',0.3,'kerneltype',ktype,'kernelpar',kpar));
setups{4} = kap(struct('p',3,'eta',.01,'eps',7e-2,'mu0',0.3,'kerneltype',ktype,'kernelpar',kpar));
setups{5} = aldkrls(struct('nu',.7,'kerneltype',ktype,'kernelpar',kpar));

numsim = 20;

%% RUN ALGORITHMS
t1 = tic;
fprintf('Fig. 2 from "Online Prediction of Time Series Data With Kernels".\n');

num_setup = length(setups);
MSE = zeros(N,num_setup);
titles = cell(num_setup,1);

for sim_ind = 1:numsim,
    fprintf('SIM %d:\n',sim_ind)
   
    % Generate the data
    [X,y,yref] = generate_richardbench(N);
    
    for setup_ind = 1:num_setup,
        kaf = setups{setup_ind};
        titles{setup_ind} = upper(class(kaf));
        
        for n=1:N
            if ~mod(n,floor(N/10)), fprintf('.'); end % progress indicator, 10 dots
            
            y_est = kaf.evaluate(X(n,:)); % test on test set
            err = yref(n) - y_est;
            MSE(n,setup_ind) = MSE(n,setup_ind) + err.^2/numsim;
            
            kaf = kaf.train(X(n,:),y(n)); % train with one input-output pair
        end
        fprintf('\n');
    end
end

toc(t1)
%% OUTPUT

% MSE smoothing by moving average for visualization
MSE_smooth = filter(1/20*ones(20,1),1,MSE);

figure
plot(10*log10(MSE_smooth));
title('Learning curves')
grid on

xlabel('iteration')
ylabel('MSE (dB)')
legend(titles)
