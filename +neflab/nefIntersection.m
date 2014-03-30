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


    inFile = [tempdir 'nefTemp.txt'];
    outFile = [tempdir 'nefOut.txt'];

    fh = fopen(inFile, 'w');
    neflab.writeMultiOFF(fh, v1, f1);
    neflab.writeMultiOFF(fh, v2, f2);
    fclose(fh);

    if ~ismac
        cmd = sprintf('env -u LD_LIBRARY_PATH NefLab intersection < %s > %s', ...
            inFile, outFile);
    else
        cmd = sprintf('NefLab intersection < %s > %s', inFile, outFile);
    end

    [status stdout] = unix(cmd);

    fh = fopen(outFile, 'r');
    [vertices faces] = neflab.readNefPolyhedron(fh);
    fclose(fh);


end