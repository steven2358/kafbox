% UNIT_TEST_PROFILER Test function for profilers of kernel adaptive
% filtering algorithms.
%
% The input arguments contains the algorithms to be tested, as separate
% strings. If no argument is provided, or a single 'all' argument, all
% algorithms from the 'lib' folder are tested.
% USAGE: unit_test_profiler('klms')
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function unit_test_profiler(varargin)

% get list of algorithms in 'lib' folder
pathstr = fileparts(which('unit_test_profiler.m'));
fdir = [pathstr '/../profiler']; % profiler folder
files = dir(fullfile(fdir,'*.m'));
[~,allfiles] = cellfun(@fileparts, {files.name}, 'UniformOutput',0);

% remove files that do not represent classes
allfiles = filter_classes(allfiles);

% get list of algorithms to run
if isempty(varargin) || (length(varargin)==1 && strcmp(varargin{1},'all'))
    algorithms = allfiles;
else
    algorithms = lower(varargin);
    for i=1:length(algorithms)
        algorithms{i} = [algorithms{i} '_profiler'];
    end
    
    fdir = ''; % local folder
    files = dir(fullfile(fdir,'*.m'));
    [~,localfiles] = cellfun(@fileparts, {files.name}, 'UniformOutput',0);
    % remove files that do not contain classes
    localfiles = filter_classes(localfiles);
    allfiles = horzcat(allfiles,localfiles);
end

% check for invalid names
for a=algorithms,
    unexisting = isempty(find(ismember(allfiles,a{1}),1));
    if unexisting
        error('Algorithm not found: %s.',a{1});
    end
end

% perform test for each specified algorithm
fprintf('\n')
for ii=1:length(algorithms)
    t1 = tic;
    algorithm = algorithms{ii};
    
    fprintf('%d. %s:\n',ii,upper(algorithm));
    fprintf('Constructing object............\n')
    try
        feval(algorithm);
    catch err
        error('Error: %s object cannot be constructed: %s',...
            upper(algorithm), err.message);
    end
    
    fprintf('Testing arguments..............\n') % constructor test 1
    opts = struct();
    try
        feval(algorithm, opts);
    catch err
        me = err.stack(1);
        error(['Error in '...
            '<a href="matlab: opentoline(which(''%s''),1)">%s</a>',...
            ' constructor: ',...
            '<a href="matlab: opentoline(which(''%s''),%d)">',...
            '%s at %d</a>\n',...
            '%s'],...
            algorithm, algorithm, me.file, me.line, me.name, me.line,...
            err.message);
    end
    
    % constructor test 2: test if bogus arguments are processed
    opts = struct('shmaloombah',5,'ker',4);
    try
        kaf = feval(algorithm, opts);
    catch err
        me = err.stack(1);
        error(['Constructor of %s object is prone ',...
            'to bogus arguments: ',...
            '<a href="matlab: opentoline(which(''%s''),%d)">',...
            '%s at %d</a>\n',...
            '%s'],...
            upper(algorithm), me.file, me.line, me.name, me.line,...
            err.message);
    end
    
    % generate some data
    c = 10;
    N = 1000;
    x = rand(N,2)*c;
    y = sin(3*x(:,1)).*cos(x(:,1)+x(:,2));
    
    fprintf('1 training step................\n')
    for i=1,
        try
            kaf = kaf.train_profiled(x(i,:),y(i));
        catch err
            me = err.stack(1);
            %             keyboard
            error(['Error in '...
                '<a href="matlab: opentoline(which(''%s''),1)">%s</a>',...
                ' training: ',...
                '<a href="matlab: opentoline(which(''%s''),%d)">',...
                '%s at %d</a>\n',...
                '%s'],...
                algorithm, algorithm, me.file, me.line, me.name, me.line,...
                err.message);
        end
    end
    
    fprintf('Profiled training....')
    flops_all = zeros(N,1);
    bytes_all = zeros(N,1);
    for i=1:N,
        if ~mod(i,floor(N/10)), fprintf('.'); end
        kaf = kaf.train_profiled(x(i,:),y(i));
        
        flops_all(i) = kaf.lastflops();
        bytes_all(i) = kaf.lastbytes();
    end
    fprintf('\n')
    
    close all;
    figure; plot(bytes_all); ylabel('bytes')
    figure; plot(flops_all); ylabel('flops')
    drawnow;
    
    t2 = toc(t1);
    fprintf('OK. %.2fs\n\n',t2);
end

if length(algorithms)>1
    fprintf('All OK.\n\n')
end


% remove files that do not represent classes
function files = filter_classes(files)
inds = [];
for i=1:length(files)
    if ~exist(files{i},'class')
        inds = [inds i]; %#ok<AGROW>
    end
end
files(inds) = [];
