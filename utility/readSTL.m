function [vertices faces] = readSTL(fname)
% [vertices faces] = readSTL(fname)

[numVerts numFaces] = scanSTL(fname);

[redundantVertices, faces] = loadData(fname, numVerts, numFaces);

[vertices, faces] = removeDuplicateVertices(redundantVertices, faces);

end

function [v f] = removeDuplicateVertices(redundantVertices, faces)

[v iUnique iInvert] = unique(redundantVertices, 'rows');

% As a result:
% v = redundantVertices(iUnique)
% redundantVertices == v(iInvert)
%
% Use to replace redundant vertex indices in the faces array.

f = iInvert(faces);

end

function [redundantVertices, faces] = loadData(fname, numVerts, numFaces)

fh = fopen(fname, 'r');

done = 0;

redundantVertices = zeros(numVerts, 3);
faces = zeros(numFaces, 3);

fNum = 1;
vNum = 1;
lineNum = 1;
while false == done
    
    line = fgetl(fh);
    
    token = strtok(line);
    
    switch token
        case 'solid'
        case 'facet'
        case 'outer'
        case 'vertex'
            [dat count] = sscanf(line, ' vertex %f %f %f');
            
            if count ~= 3
                poop(lineNum, line, count);
            end
            
            redundantVertices(vNum,:) = dat;
            vNum = vNum + 1;
            
        case 'endloop'
        case 'endfacet'
            faces(fNum,:) = vNum-3:vNum-1;
            fNum = fNum + 1;
            
        case 'endsolid'
            done = 1;
        otherwise
            poop(lineNum, line);
    end
    
    lineNum = lineNum + 1;
end

if fNum ~= numFaces + 1
    warning('Found %i faces, expecting %i', fNum-1, numFaces);
end

if vNum ~= numVerts + 1
    warning('Found %i vertices, expecting %i', vNum-1, numVerts);
end



fclose(fh);

end


function poop(lineNum, line, count)

if nargin == 2
    error('Cannot understand line %i: %s', lineNum, line);
elseif nargin == 3
    error('Only parsed %i tokens of line %i: %s', count, lineNum, line);
end

end


function [numVerts numFaces] = scanSTL(filename)

fh = fopen(filename, 'r');

numFaces = 0;
numVerts = 0;

done = false;
while false == done
    line = fgetl(fh);
    token = strtok(line);
    
    if strcmpi(token, 'facet')
        numFaces = numFaces + 1;
    elseif strcmpi(token, 'vertex')
        numVerts = numVerts + 1;
    elseif strcmpi(token, 'endsolid')
        done = true;
    end
    
end

fclose(fh);

end
