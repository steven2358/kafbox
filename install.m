% KAFBOX installation file. Adds local folders to path.

disp('Adding KAFBOX folders to Matlab path...')

addpath(genpath([pwd '/data']));
addpath(genpath([pwd '/demo']));
addpath(genpath([pwd '/lib']));
addpath(genpath([pwd '/lib/util']));

disp('Saving Matlab path...')

savepath;