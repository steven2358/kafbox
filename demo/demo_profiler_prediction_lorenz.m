% Demo: demonstration of the kernel adaptive filter algorithm profiler. 
% Compares the cost vs prediction error tradeoffs and convergence speeds
% for several algorithms on the Lorenz data set.

clear all
close all

%% PARAMETERS

output_dir_default = 'profiler_output';

% data and algorithm setup
data.name = 'Lorenz';
data.file = 'lorenz.dat';
data.N = 10000; % number of data points
data.embedding = 6; % time embedding
data.numsim = 5; % 10 minutes per simulation on an Intel Pentium Core2 Duo
data.offset = 50; % apply offset per simulation
data.error_measure = 'MSE';

i=0; % initialize setups

%% ALD-KRLS
i=i+1;
algorithms{i}.name = 'ALD-KRLS';
algorithms{i}.class = 'aldkrls';
algorithms{i}.figstyle = struct('color',[.75  0 .75],'marker','^');
algorithms{i}.options = struct('sweep_par','nu','sweep_val',[1E-4 2E-4 1E-3 .01 .05 .1],...
    'kerneltype','gauss','kernelpar',32);

%% Q-KLMS
i=i+1;
algorithms{i}.name = 'Q-KLMS';
algorithms{i}.class = 'qklms';
algorithms{i}.figstyle = struct('color',[1  0  0],'marker','o');
algorithms{i}.options = struct('eta',0.5,'sweep_par','epsu','sweep_val',[1 2 5 10 12 15 18],...
    'kerneltype','gauss','kernelpar',32);

%% SW-KRLS
i=i+1;
algorithms{i}.name = 'SW-KRLS';
algorithms{i}.class = 'swkrls';
algorithms{i}.figstyle = struct('color',[0  .75 .75],'marker','s');
algorithms{i}.options = struct('c',1E-6,'sweep_par','M','sweep_val',[1 2 3 4 5 10 20 50 100 200],...
    'kerneltype','gauss','kernelpar',32);

%% FB-KRLS
i=i+1;
algorithms{i}.name = 'FB-KRLS';
algorithms{i}.class = 'fbkrls';
algorithms{i}.figstyle = struct('color',[0  .5  0],'marker','d');
algorithms{i}.options = struct('lambda',1E-6,'sweep_par','M','sweep_val',[10 15 20 25 30 50 100],...
    'kerneltype','gauss','kernelpar',32);

%% KRLS-T
i=i+1;
algorithms{i}.name = 'KRLS-T';
algorithms{i}.class = 'krlst';
algorithms{i}.figstyle = struct('color',[0  0  1],'marker','+');
algorithms{i}.options = struct('sn2',1E-6,'lambda',1,'sweep_par','M','sweep_val',[10 15 17 20 24 30 50 100 200],...
    'kerneltype','gauss','kernelpar',32);

%% PROGRAM

fprintf('Running profiler for %d algorithms on %s data.\n',i,data.name);
output_dir = input('Folder for storing results (ENTER = subfolder here): ','s');
if isempty(output_dir),
    output_dir = output_dir_default;
    fprintf('Using default folder "%s" for storing results.\n', output_dir_default);
end

t1 = tic;
[data,algorithms,results] = kafbox_profiler(data,algorithms,output_dir);
t2 = toc(t1);

fprintf('Elapsed time: %d seconds\n',ceil(t2));

%% OUTPUT

mse_curves = kafbox_profiler_msecurves(results);

resinds = [1,1;2,1;3,10;4,7;5,9]; % result indices
[f0,h0] = kafbox_profiler_plotconvergence(algorithms,mse_curves,resinds);

[f1,h1] = kafbox_profiler_plotresults(algorithms,mse_curves,results,{'bytes','flops'});
object_handle = get(h1);
legend(object_handle.String,'Location','NW'); % move legend

[f2,h2] = kafbox_profiler_plotresults(algorithms,mse_curves,results,{'ssmse','flops'});

[f3,h3] = kafbox_profiler_plotresults(algorithms,mse_curves,results,{'ssmse','bytes'});

[f4,h4] = kafbox_profiler_plotresults(algorithms,mse_curves,results,{'ssmse','ssplus'},1);
object_handle = get(h4);
legend(object_handle.String,'Location','SW'); % move legend

[f5,h5] = kafbox_profiler_plotresults(algorithms,mse_curves,results,{'ssmse','timeto'},-20);
object_handle = get(h5);
legend(object_handle.String,'Location','NW'); % move legend

[f6,h6] = kafbox_profiler_plotresults(algorithms,mse_curves,results,{'ssmse','timetopct'},90);
object_handle = get(h6);
legend(object_handle.String,'Location','SW'); % move legend
