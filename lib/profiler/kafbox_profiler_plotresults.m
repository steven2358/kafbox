% Produces a plot showing the performance trade-off for two chosen
% performance measures for different algorithms. Includes a data cursor
% showing extra data.
%
% Input:
% - algorithms: cell containing algorithm setups
% - mse_curves: cell containing MSE curves (possibly averaged out)
% - results: cell containing complete results
% - pmeasures: 2x1 cell containing two performance measures
% - opt: optional parameters / options for performance measures
%
% Output:
% - f: figure axis handle
% - h: legend axis handle
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function [f,h] = kafbox_profiler_plotresults(algorithms,mse_curves,...
    results,pmeasures,opt)

logaxis = [0 0];
labels = cell(2);
for i=1:2
    switch pmeasures{i}
        case 'bytes'
            labels{i} = 'max bytes';
            logaxis(i) = 1;
        case 'flops'
            labels{i} = 'max flops';
            logaxis(i) = 1;
        case 'ssmse'
            labels{i} = 'steady-state MSE';
        case 'ssplus'
            labels{i} = sprintf('time to ss+%d dB',opt);
            logaxis(i) = 1;
        case 'timeto'
            labels{i} = sprintf('time to %d dB',opt);
            logaxis(i) = 1;
        case 'timetopct' % time to fall to percentage of steady state
            labels{i} = sprintf('time to %d%% of steady-state',opt);
            logaxis(i) = 1;
        otherwise
            error('undefined case')
    end
end

f = figure; hold all %#ok<NASGU>
set(gcf,'Position',[200, 200, 500 300])

starts = zeros(length(algorithms),2);
titles = cell(length(algorithms),1);
for algo_ind = 1:length(algorithms)
    algo = algorithms{algo_ind};
    
    algo_results = results{algo_ind};
    xdata = zeros(length(algo_results),1);
    ydata = zeros(length(algo_results),1);
    
    num_algo = length(algo_results);
    for config_ind = 1:num_algo
        config_results = algo_results{config_ind};
        
        numsim = length(config_results);
        flops = 0;
        bytes = 0;
        for sim_ind = 1:numsim,
            sim_results = config_results{sim_ind};
            bytes = max(bytes,sim_results.bytes);
            flops = max(flops,sim_results.flops);
        end
        
        mse_curve = mse_curves{algo_ind}{config_ind};
        ss = kafbox_profiler_convergence_analysis(mse_curve);
        r.flops = flops;
        r.bytes = bytes;
        r.ssmse = ss;
        % keyboard
        clear target_mse
        for i=1:2,
            switch pmeasures{i}
                case 'ssplus'
                    target_mse = ss + opt;
                case 'timeto'
                    target_mse = opt;
                case 'timetopct'
                    target_mse = ss * opt/100;
            end
            if exist('target_mse','var')
                [ss,tt] = kafbox_profiler_convergence_analysis(...
                    mse_curve,target_mse);
                r.(pmeasures{i}) = tt;
            end
        end
        
        xdata(config_ind) = r.(pmeasures{1});
        ydata(config_ind) = r.(pmeasures{2});
    end
    
    %starts(algo_ind,:) = [xdata(1),ydata(1)];
    
    titles{algo_ind} = algo.name;
    
    % tooltip data
    uservalues = algo.options.sweep_val;
    userstrings.method = algo.name;
    userstrings.param = algo.options.sweep_par;
    userstrings.x = xdata;
    userstrings.y = ydata;
    userstrings.xaxis = labels{1};
    userstrings.yaxis = labels{2};
    
    
    xxi = xdata; if logaxis(1), xxi = log10(xdata); end
    yyi = ydata; if logaxis(2), yyi = log10(ydata); end
    starts(algo_ind,:) = [xxi(1),yyi(1)];
    hLine = line(xxi, yyi,...
        'Marker', algo.figstyle.marker,...
        'MarkerSize', 8,...
        'Color', algo.figstyle.color,...
        'UserData', {uservalues, userstrings},...
        'LineStyle','-',...
        'LineWidth',2); %#ok<NASGU>
    
    cursorMode = datacursormode(gcf);
    set(cursorMode, 'enable','on', 'UpdateFcn',...
        @setDataTipTxt, 'NewDataCursorOnClick',false);
    
end
h = legend(titles);
grid on; box on
xlabel(labels{1});
ylabel(labels{2});

if logaxis(1),
    xtl = get(gca,'XTickLabel');
    set(gca,'XTickLabel',str2log10(xtl));
end
if logaxis(2),
    ytl = get(gca,'YTickLabel');
    set(gca,'YTickLabel',str2log10(ytl));
end

% Expand figure axes
f = gcf;
style = hgexport('factorystyle'); style.Bounds = 'tight';
hgexport(f,'-clipboard',style,'applystyle', true);
drawnow;

% plot start points
for algo_ind=1:length(algorithms)
    xi = starts(algo_ind,1);
    yi = starts(algo_ind,2);
    yyi = yi;
    hLine = line(xi(1), yyi(1),...
        'Marker', 'o', 'MarkerSize', 6,...
        'Color', 'black', 'MarkerFaceColor','black'); %#ok<NASGU>
end


function output_txt = setDataTipTxt(~,event_obj)
% ~            Currently not used (empty)
% event_obj    Object containing event data structure
% output_txt   Data cursor text (string or cell array of strings)

p = event_obj.Position;
t = event_obj.Target;
ud = get(t,'UserData');

% get index in data series by looking up closest point
p_all = [get(t,'XData')', get(t,'YData')'];
N = size(p_all,1);
d = p_all-repmat(p,N,1);
d2 = sum(d.^2,2);
[mm,ii] = min(d2); %#ok<ASGLU>

uservalue = ud{1}(ii);
userstrings = ud{2};

% fill text tooltip
output_txt{1} = userstrings.method;
output_txt{2} = print_formatted_value(userstrings.param,uservalue);
output_txt{3} = print_formatted_value(userstrings.xaxis,userstrings.x(ii));
output_txt{4} = print_formatted_value(userstrings.yaxis,userstrings.y(ii));


function str = print_formatted_value(s,v)
if (v-round(v))^2<1E-8,
    str = sprintf('%s = %d',s,v);
else
    str = sprintf('%s = %.2f',s,v);
end

% add log10 to string numbers
function b = str2log10(a)
b = cell(length(a),1);
for i=1:length(a),
    b{i} = ['1E' strtrim(a(i,:))];
end
