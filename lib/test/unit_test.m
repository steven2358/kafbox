% UNIT_TEST Test function for kernel adaptive filtering algorithms.
%
% The input arguments contains the algorithms to be tested, as separate
% strings. If no argument is porovided, or a single 'all' argument, all
% algorithms are tested.
% USAGE: unit_test('klms','krlst')
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% http://sourceforge.net/projects/kafbox/

function unit_test(varargin)

% get list of all algorithms
files = dir(fullfile('../','*.m'));
[~,files] = cellfun(@fileparts, {files.name}, 'UniformOutput',false);

% convert into list of algorithms
if isempty(varargin) || (length(varargin)==1 && strcmp(varargin{1},'all'))
    algorithms = files;
else
    algorithms = lower(varargin);
end

% check for invalid names
for a=algorithms,
    unexisting = isempty(find(ismember(files,a{1}),1));
    if unexisting
        error('Wrong algorithm acronym: %s',a{1});
    end
end

fprintf('\n')
for ii=1:length(algorithms)
    algorithm = algorithms{ii};
    
    fprintf('%d. %s:\n',ii,upper(algorithm));
    fprintf('Constructing object...\n')
    kaf = feval(algorithm);

    fprintf('Testing training phase')
    c = 10;
    N = 10000;
    x = rand(N,2)*c;
    y = sin(3*x(:,1)).*cos(x(:,1)+x(:,2));
    for i=1:N-1,
        if ~mod(i,floor(N/10)), fprintf('.'); end
        kaf = kaf.train(x(i,:),y(i));
        % y_test = kaf.evaluate(x(i+1,:));
    end
    fprintf('\n')

    fprintf('Testing evaluation phase...\n')
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
    
fprintf('All OK.\n\n')