function [vertices faces] = nefUnion(v1, f1, v2, f2)
% nefUnion    Calculate union of two polyhedra
%
% [vertices faces] = nefUnion(v1, f1, v2, f2)
%

% Write args to file
% Call NefLab
% Get polyhedron back out

inFile = [tempdir 'nefTemp.txt'];
outFile = [tempdir 'nefOut.txt'];

fh = fopen(fname, 'w');
neflab.writeMultiOFF(fh, v1, f1);
neflab.writeMultiOFF(fh, v2, f2);
fclose(fh);

[status, stdout] = ...
    unix(sprintf('NefLab union < %s > %s', inFile, outFile));

fh = fopen(outFile, 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

