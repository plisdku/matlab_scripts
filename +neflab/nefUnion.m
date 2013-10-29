function [vertices faces] = nefUnion(v1, f1, v2, f2)
% nefUnion    Calculate union of two polyhedra
%
% [vertices faces] = nefUnion(v1, f1, v2, f2)
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
    unix('env -u LD_LIBRARY_PATH NefLab union < nefTemp.txt > nefOut.txt');

fh = fopen('nefOut.txt', 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

