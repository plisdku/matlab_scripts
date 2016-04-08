function tris = triangulate(vertices, varargin)

numContours = numel(varargin);

if numContours == 1 && numel(varargin{1}) == 3
    % It's already a triangle, just accept it as-is.
    tris = reshape(varargin{1}, 1, 3);
    return
else
    [local2global, localContours, contourLengths] = ...
        localIndicesRobust(varargin{:});
    localVertices = vertices(local2global,:);
    
    % Beginning and ending indices of contours, e.g.
    %   cStarts = [1 5 8]
    %   cEnds = [4 7 15]
    % means the outer contour is localVertices(1:4,:)
    % and the inner contours run from 5 to 7 and from 8 to 15.
    
    cStarts = 1 + cumsum([0 contourLengths(1:end-1)]);
    cEnds = cStarts + contourLengths - 1;
    
    flatVertices = flatten(localVertices, ...
        neflab.facetNormal(localVertices, ...
            localContours(cStarts(1):cEnds(1))));
    
    % Explicitly list all the edges of the polygon, for Delaunay
    % triangulation.  These edges are triangulation constraints.
    col = @(A) reshape(A, [], 1);
    constraints = [col(localContours), col(localContours)];
    
    for cc = 1:numel(contourLengths)
        firstCol = localContours(cStarts(cc):cEnds(cc),1);
        secondCol = firstCol([2:end,1]);
        constraints(cStarts(cc):cEnds(cc),2) = secondCol;
    end
    
    if size(uniqueRows(flatVertices), 1) ~= size(flatVertices,1)
        disp('Duplicate vertices found')
    end
    %dt = delaunayTriangulation(flatVertices, constraints);
    dt = DelaunayTri(flatVertices, constraints);
    inside = inOutStatus(dt);
    
    localTris = dt.Triangulation(inside,:);
    tris = local2global(localTris);
    
    if size(tris,2) ~= 3
        tris = tris';
    end
    
    %{
    % It's not a triangle.  It has more vertices and might have holes.
    [localContours, local2global, uniqueGlobal] = localIndices(varargin{:});
    localVertices = vertices(uniqueGlobal,:);

    outerContour = localContours{1};

    nv = neflab.facetNormal(localVertices, outerContour);
    flatVertices = flatten(localVertices, nv);

    colVec = @(A) reshape(A, [], 1);

    constraints = [colVec(outerContour(1:end)),...
        colVec(outerContour([2:end,1]))];

    for innerLoop = 2:length(localContours)

        loop = localContours{innerLoop};

        constraints = [constraints; ...
            [colVec(loop(1:end)), colVec(loop([2:end,1]))] ];
    end
    
    dt = DelaunayTri(flatVertices, constraints);
    inside = inOutStatus(dt);

    localTris = dt.Triangulation(inside,:);

    tris = local2global(localTris);

    if size(tris,2) ~= 3
        tris = tris';
    end
    %}
end




end

%{
% Convert global vertex indices to local facet indices.
% localContours will take on values from 1 to numVerts (in facet).
% local2global gives back the global indices.
function [localContours, local2global, uniqueGlobal] = localIndices(varargin)

numContours = numel(varargin);
localContours = varargin;

numVerts = sum(cellfun(@numel, varargin));
local2global = cat(1, varargin{:});
[uniqueGlobal, ia, localContour] = unique(local2global, 'stable');

iVert = 0;
for cc = 1:numContours
    contourLength = numel(varargin{cc});
    localContours{cc} = iVert + (1:contourLength)';
    iVert = iVert + contourLength;
end

assert(numel(uniqueGlobal) == numel(unique(local2global)));

end
%}

function [local2global, localContours, contourLengths] = ...
    localIndicesRobust(varargin)

numContours = numel(varargin);
contourLengths = cellfun(@numel, varargin);

%cStarts = 1 + cumsum([0 contourLengths(1:end-1)]);
%dcEnds = cStarts + contourLengths - 1;

[local2global, ~, localContours] = unique(cat(1,varargin{:}), 'stable');

end


% The old way that could handle repeated vertices.  I want no part of that
% silliness.  It was slow.  :-)
%{
% concatenate all global vertices
globalVerts = [];
for cc = 1:numContours
    globalVerts = [globalVerts(:); varargin{cc}];
end

% Build the maps!  Thank you, Matlab, for making this so easy!!
% uniqueGlobal   the global vertex indices used in this facet
% local2global   an array that turns local indices (1 to Nv) into globals
% localIndices   a number in 1:Nv for each global index in varargin
[uniqueGlobal, ia, localIndices] = unique(globalVerts);
[~, uniqueLocal] = sort(uniqueGlobal);
local2global = globalVerts(ia);

localVerts = uniqueLocal(localIndices);

% Create the local contours now.  It's a matter of replacing indices.
localContours = varargin;
iVert = 1;

for cc = 1:numContours
    contourLength = numel(localContours{cc});
    localContours{cc} = localIndices(iVert:(iVert + contourLength - 1));
    iVert = iVert + contourLength;
end
%}


% I need to create a right-handed coordinate system in which normalVector
% is "up", then project vIn onto this perpspace of normalVector.
%
% fast method from Jeppe Revall Frisvad.  Yay Denmark!!  :-)
function vOut = flatten(vIn, normalVector)

% normalVector should be unit length already!

    if normalVector(3) < -0.9999999
        b1 = [0; -1; 0];
        b2 = [-1; 0; 0];
    else
        a = 1.0/(1.0 + normalVector(3));
        b = -normalVector(1) * normalVector(2) * a;
        b1 = [1.0 - normalVector(1)^2*a; b; -normalVector(1)];
        b2 = [b; 1.0 - normalVector(2)^2*a; -normalVector(2)];
    end
    
    perpSpace = [b1 b2];
    % make sure it's a right-handed coordinate system, and has
    % determinant 1-ish.
    assert(abs(det([perpSpace, reshape(normalVector, [], 1)]) - 1.0) < 1e-4);
    
    % project everything!
    vOut = vIn * perpSpace;
end


% This is the old, grossly slow method.
%{
colVec = @(A) reshape(A, [], 1);
rowVec = @(A) reshape(A, 1, []);

perpSpace = null(rowVec(normalVector));

if det([perpSpace, colVec(normalVector)]) < 0
    perpSpace = fliplr(perpSpace);
end

% Now perpSpace is two column vectors of right-handed orientation with
% respect to normalVector.  Thank you, Matlab!!
%
% Perform the projection.

%projectionMatrix = perpSpace';

vOut = vIn * perpSpace;
%}



