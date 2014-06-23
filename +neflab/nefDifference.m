function [vertices, faces] = nefDifference(v1, f1, v2, f2)
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

inFile = [tempdir sprintf('nefTemp%1.4f.txt', now)];
outFile = [tempdir sprintf('nefOut%1.4f.txt', now)];

fh = fopen(inFile, 'w');
neflab.writeMultiOFF(fh, v1, f1);
neflab.writeMultiOFF(fh, v2, f2);
fclose(fh);

if ~ismac
    cmd = sprintf('env -u LD_LIBRARY_PATH NefLab difference < %s > %s', ...
        inFile, outFile);
else
    cmd = sprintf('NefLab difference < %s > %s', inFile, outFile);
end

[status, stdout] = unix(cmd);

if status
    keyboard;
end

fh = fopen(outFile, 'r');
[vertices faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

