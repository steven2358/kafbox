function [f,h] = kafbox_profiler_plotconvergence(algorithms,...
    mse_curves,resinds)

figure; hold all
set(gcf,'Position',[200, 200, 500 300])

titles = cell(length(algorithms),1);

for i=1:size(resinds,1)
    algo = algorithms{resinds(i,1)};
    
    ls = '-'; % line style
    if isfield(algo.figstyle,'line')
        ls = algo.figstyle.line;
    end
    
    curve = mse_curves{resinds(i,1)}{resinds(i,2)};
    xs= ~isnan(curve);
    inds = 1:length(mse_curves{resinds(i,1)}{resinds(i,2)});
    plot(inds(xs),10*log10(curve(xs)),'color',algo.figstyle.color,...
        'LineWidth',1,'LineStyle',ls)
    
    titles{i} = algo.name;
end

h = legend(titles);
grid on; box on
xlabel('iteration');
ylabel('(N)MSE');

% Expand figure axes
f = gcf;
style = hgexport('factorystyle'); style.Bounds = 'tight';
hgexport(f,'-clipboard',style,'applystyle', true);
drawnow;
