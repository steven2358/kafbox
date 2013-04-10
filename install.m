% KAFBOX installation file. Adds local folders to path.

disp('Adding KAFBOX folders to MATLAB path...')

addpath(genpath([pwd '/demo']));
addpath(genpath([pwd '/lib']));

savepath;