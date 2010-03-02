% materialspath.m
function materialspath

sep = filesep;
thisfile = mfilename('fullpath');

[pathstr, name] = fileparts(thisfile);

fdtd_parent = fileparts(pathstr);

pathlist{1} = pathstr;

for nn = 1:length(pathlist)
    addpath(pathlist{nn});
end

disp('The Matlab path has been modified to include this script distribution.');
disp('In order to have these scripts permanently available in your path,');

if verLessThan('matlab', '7')
    theTool = '<a href="matlab:pathtool">pathtool</a>';
else
    theTool = '<a href="matlab:savepath">savepath</a>';
end

disp(sprintf('run %s from the command line or modify startup.m with the lines\n', ...
    theTool));

for nn = 1:length(pathlist)
    disp(sprintf('addpath %s', pathlist{nn}));
end

