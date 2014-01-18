% Proces MSE results into MSE curves: average out and/or limit to a
% subset of temporal indices.
%
% Input:
% - results: a configresults cell produced by kafbox_profiler.
% - inds: array of indices for which to calculate output
%
% Output:
% - MSE_avg_setups: cell containing averaged out MSE curves, structure
% corresponds to "results" cell structure.
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function MSE_avg_setups = kafbox_profiler_msecurves(results,inds)

num_setups = length(results);
MSE_avg_setups = cell(num_setups,1);

for setup_ind = 1:num_setups,
    setup_results = results{setup_ind};
    
    num_configs = length(setup_results);
    MSE_avg_configs = cell(num_configs,1);

    for config_ind = 1:num_configs,
        config_results = setup_results{config_ind};
        
        N = length(config_results{1}.MSE);
        MSE = zeros(N,1);
        
        num_sim = length(config_results);
        for sim_ind = 1:num_sim,
            simresults = config_results{sim_ind};
            
            if isfield(simresults,'NMSE')
                simresults.MSE = simresults.NMSE;
            end
            if nargin<2
                inds = 1:length(simresults.MSE);
            end

            mm = min(length(simresults.MSE(inds)),N);
            MSE = MSE(1:mm) + simresults.MSE(inds(1:mm))/num_sim;
        end
        MSE_avg_configs{config_ind} = MSE;
    end
    MSE_avg_setups{setup_ind} = MSE_avg_configs;
end
