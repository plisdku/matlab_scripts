function [indicesNew, indicesOld] = paternity(v1, f1, v2, f2)
% paternity    Calculate which faces of polyhedron 2 came from which faces
% of polyhedron 1
%
% [indicesNew, indicesOld] = paternity(v1, f1, v2, f2)
%
% 

% Write args to file
% Call NefLab
% Get polyhedron back out

if neflab.disjointHulls(v1, v2)
    vertices = v1;
    faces = f1;
    return;
end

inFile = [tempdir sprintf('nefTemp%1.7f.txt', now)];
outFile = [tempdir sprintf('nefOut%1.7f.txt', now)];

fh = fopen(inFile, 'w');
assert(fh ~= -1);
neflab.writeMultiOFF(fh, v1, f1);
neflab.writeMultiOFF(fh, v2, f2);
fclose(fh);

if ~ismac
    cmd = sprintf('env -u LD_LIBRARY_PATH NefLab inherit < %s > %s', ...
        inFile, outFile);
else
    cmd = sprintf('NefLab inherit < %s > %s', inFile, outFile);
end

[status, stdout] = unix(cmd);

if status
    keyboard;
end

[fh, message] = fopen(outFile, 'r');
assert(fh ~= -1);
[indicesNew, indicesOld] = readPaternity(fh);
%[vertices, faces] = neflab.readNefPolyhedron(fh);
fclose(fh);

end

function [indicesNew, indicesOld] = readPaternity(fh)

indicesNew = [];
indicesOld = [];

assert(fh ~= -1);

iNew = 1;
while 1
    tline = fgetl(fh);
    if ~ischar(tline)
        break
    end
    
    [~, ~, ~, nextindex] = sscanf(tline, '%d:', 1);
    [indices, count, errmsg, nextindex] = sscanf(tline(nextindex:end), '%d', [1,inf]); 
    
    indicesNew = [indicesNew; iNew*ones(length(indices),1)];
    indicesOld = [indicesOld; reshape(indices+1, [], 1)]; % convert C to Matlab indices
    
    iNew = iNew + 1;
    
end

end