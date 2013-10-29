function [vertices faces] = nefDifference(v1, f1, v2, f2)
% nefDifference    Calculate difference of two polyhedra
%
% [vertices faces] = nefDifference(v1, f1, v2, f2)
%

% Write args to file
% Call NefLab
% Get polyhedron back out

if neflab.disjointHulls(v1, v2)
    vertices = v1;
    faces = f1;
    return;
end

fname = 'nefTemp.txt';

fh = fopen(fname, 'w');
neflab.writeMultiOFF(fh, v1, f1);
neflab.writeMultiOFF(fh, v2, f2);
fclose(fh);

[status, stdout] = ...
    unix('env -u LD_LIBRARY_PATH NefLab difference < nefTemp.txt > nefOut.txt');

fh = fopen('nefOut.txt', 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

