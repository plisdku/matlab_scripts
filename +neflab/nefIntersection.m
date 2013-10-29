function [vertices faces] = nefIntersection(v1, f1, v2, f2)
% nefIntersection    Calculate intersection of two polyhedra
%
% [vertices faces] = nefIntersection(v1, f1, v2, f2)
%

% Write args to file
% Call NefLab
% Get polyhedron back out

if neflab.disjointHulls(v1, v2)
    vertices = [];
    faces = [];
    return;
end

fname = 'nefTemp.txt';

fh = fopen(fname, 'w');
neflab.writeMultiOFF(fh, v1, f1);
neflab.writeMultiOFF(fh, v2, f2);
fclose(fh);

[status, stdout] = ...
    unix('env -u LD_LIBRARY_PATH NefLab intersection < nefTemp.txt > nefOut.txt');

fh = fopen('nefOut.txt', 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

