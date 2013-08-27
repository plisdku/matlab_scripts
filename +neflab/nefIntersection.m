function [vertices faces] = nefIntersection(v1, f1, v2, f2)
% nefIntersection    Calculate intersection of two polyhedra
%
% [vertices faces] = nefIntersection(v1, f1, v2, f2)
%

% Write args to file
% Call NefLab
% Get polyhedron back out

fname = 'nefTemp.txt';

fh = fopen(fname, 'w');
neflab.writePolyhedron(v1, f1, fh);
neflab.writePolyhedron(v2, f2, fh);
fclose(fh);

[status, stdout] = ...
    unix('./NefLab intersection < nefTemp.txt > nefOut.txt');

fh = fopen('nefOut.txt', 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

