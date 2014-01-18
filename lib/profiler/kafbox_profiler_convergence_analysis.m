% Convergence analysis. Returns steady state MSE, numbers of iterations to
% reach certain MSE values and MSEs reached after certain numbers of
% iterations.
%
% Input:
% - results: either a configresults cell or a simresults structure produced
%   by kafbox_profiler.
% - target_mse: target mse value
%
% Output:
% - ss: steady-state MSE level over last 1000 samples
% - mse_marks: structure containing
%       - labels: array indicating target mse values
%       - iterations: array with numbers of iterations to reach these mses
% - iter_marks: structure containing
%       - labels: array indicating iteration numbers values
%       - MSEs: array with mses reached after these iterations
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function [ss,it_reached,mse_reached] = ...
    kafbox_profiler_convergence_analysis(MSE_curve, target_mse, target_it)

xs = find(~isnan(MSE_curve));

% calculate steady state as MSE over last 1000 points
ss = 10*log10(mean(MSE_curve(xs(xs>length(MSE_curve)-1000))));

if nargin > 1
    % check when error_measure reaches certain MSE
    xsi = find(10*log10(MSE_curve(xs)) < target_mse);
    if numel(xsi)>0
        if xsi(1)==1
            index = xsi(1);
        else
            % linear interpolation to find value
            sti = xsi(1)-1;
            ndi = xsi(1);
            stv = 10*log10(MSE_curve(xs(sti)));
            ndv = 10*log10(MSE_curve(xs(ndi)));
            index = xs(sti)+(xs(ndi) - xs(sti))*(stv-target_mse)/(stv-ndv);
            %     plot(index,target_MSE,'ro')
            %     plot([0 index],[target_MSE target_MSE],'r','LineWidth',2)
            %     ylim = get(gca,'YLim');
            %     plot([index index],[ylim(1) target_MSE],'r')
        end
    else
        % fprintf('No convergence to %d dB.\n',target_MSE);
        index = nan;
    end
    it_reached = index;
end

if nargin > 2
    % check what error is reached after certain number of iterations
    xsi = find((xs) >= target_it);
    if xsi(1)==1
        value = 10*log10(MSE_curve(xsi(1)));
    else
        % linear interpolation to find value
        sti = xsi(1)-1;
        ndi = xsi(1);
        stv = 10*log10(MSE_curve(xs(sti)));
        ndv = 10*log10(MSE_curve(xs(ndi)));
        value = stv + (ndv - stv)*(xs(sti)-target_iter)/(xs(sti)-xs(ndi));
        % plot(target_iter,value,'ko')
        % plot([0 target_iter],[value value],'k')
        % ylim = get(gca,'YLim');
        % plot([target_iter target_iter],[ylim(1) value],'k','LineWidth',2)
    end
    mse_reached = value;
end
