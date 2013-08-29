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
neflab.writeMultiOFF(fh, v1, f1);
neflab.writeMultiOFF(fh, v2, f2);
%neflab.writePolyhedron(v1, f1, fh);
%neflab.writePolyhedron(v2, f2, fh);
fclose(fh);

[status, stdout] = ...
    unix('NefLab union < nefTemp.txt > nefOut.txt');

fh = fopen('nefOut.txt', 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

