% UNIT_TEST Test function for kernel adaptive filtering algorithms.
%
% The input arguments contains the algorithms to be tested, as separate
% strings. If no argument is provided, or a single 'all' argument, all
% algorithms are tested.
% USAGE: unit_test('klms','krlst')
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function unit_test(varargin)

% get list of algorithms to run
if isempty(varargin) || (length(varargin)==1 && strcmp(varargin{1},'all'))
    pathstr = fileparts(which('unit_test.m'));
    fdir = [pathstr '/..']; % lib folder
    files = dir(fullfile(fdir,'*.m'));
    [~,files] = cellfun(@fileparts, {files.name}, 'UniformOutput',false);
    algorithms = files;
else
    fdir = ''; % local folder
    files = dir(fullfile(fdir,'*.m'));
    [~,files] = cellfun(@fileparts, {files.name}, 'UniformOutput',false);
    algorithms = lower(varargin);
end

% check for invalid names
for a=algorithms,
    unexisting = isempty(find(ismember(files,a{1}),1));
    if unexisting
        error('Algorithm not found: %s.',a{1});
    end
end


% perform test for each specified algorithm
fprintf('\n')
for ii=1:length(algorithms)
    algorithm = algorithms{ii};
    
    fprintf('%d. %s:\n',ii,upper(algorithm));
    fprintf('Constructing object............\n')
    try
        kaf = feval(algorithm);
    catch err
        error('Error: %s object cannot be constructed: %s',...
            upper(algorithm), err.message);
    end
    
    fprintf('Testing arguments..............\n') % bogus arguments
    opts = struct();
    try
        kaf = feval(algorithm, opts);
    catch err
        me = err.stack(1);
        error(['Error in %s constructor: ',...
            '<a href="matlab: opentoline(which(''%s''),%d)">',...
            '%s at %d</a>\n',...
            '%s'],...
            upper(algorithm), me.file, me.line, me.name, me.line,...
            err.message);
    end
    
    % test if bogus args are processed
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
    N = 10000;
    x = rand(N,2)*c;
    y = sin(3*x(:,1)).*cos(x(:,1)+x(:,2));
    
    fprintf('Testing 1 training step........\n')
    for i=1,
        try
            kaf = kaf.train(x(i,:),y(i));
        catch err
            me = err.stack(1);
            error(['Error in %s training: ',...
                '<a href="matlab: opentoline(which(''%s''),%d)">',...
                '%s at %d</a>\n',...
                '%s'],...
                upper(algorithm), me.file, me.line, me.name, me.line,...
                err.message);
        end
    end
    
    fprintf('Testing long training')
    for i=1:N,
        if ~mod(i,floor(N/10)), fprintf('.'); end
        kaf = kaf.train(x(i,:),y(i));
        % y_test = kaf.evaluate(x(i+1,:));
    end
    fprintf('\n')
    
    fprintf('Testing evaluation.............\n')
    [x1,x2] = meshgrid(0:.2:c, 0:.2:c);
    yt = kaf.evaluate([x1(:) x2(:)]);
    
    z = reshape(yt,size(x1,1),size(x2,1));
    close; figure;
    surf(x1,x2,z);
    title(sprintf('%s',upper(algorithm)));
    colormap('spring'); view(20,70);
    drawnow;
    
    fprintf('OK.\n\n');
end

if length(algorithms)>1
    fprintf('All OK.\n\n')
end
