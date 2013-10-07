% Installation file. Adds local folders to path.

fprintf('Adding KAFBOX folders to Matlab path... ')

addpath(genpath([pwd '/data']));
addpath(genpath([pwd '/demo']));
addpath(genpath([pwd '/lib']));

fprintf('done.\n')
disp('Type "savepath" if you wish to store changes.')
% savepath;
