function tris = triangulate(vertices, varargin)

[localContours local2global uniqueGlobal] = localIndices(varargin{:});
localVertices = vertices(uniqueGlobal,:);

outerContour = localContours{1};

nv = neflab.facetNormal(localVertices, outerContour);

flatVertices = flatten(localVertices, nv);
%{
figure(10); clf
subplot(211)
plot3(localVertices(:,1), localVertices(:,2), localVertices(:,3), 'bo-');
subplot(212)
plot(flatVertices(:,1), flatVertices(:,2), 'bo-');
pause
%}

%[~, mainAxis] = max(abs(nv));
%flatVertices = localVertices(:,[1:mainAxis-1, mainAxis+1:3]);

colVec = @(A) reshape(A, [], 1);

constraints = [colVec(outerContour(1:end)),...
    colVec(outerContour([2:end,1]))];

for innerLoop = 2:length(localContours)
    
    loop = localContours{innerLoop};
    
    constraints = [constraints; ...
        [colVec(loop(1:end)), colVec(loop([2:end,1]))] ];
        %[colVec(loop(1:end)), colVec(loop(1:end-1))] ];
%        [colVec(loop(1:end-1)), colVec(loop(2:end))] ];
end

dt = DelaunayTri(flatVertices, constraints);
inside = inOutStatus(dt);

localTris = dt.Triangulation(inside,:);

%if numel(localContours) > 1
%    fprintf('Complicated.\n');
%end

% Go back to global vertex indices
%disp(local2global)
%pause
tris = local2global(localTris);



% Convert global vertex indices to local facet indices.
% localContours will take on values from 1 to numVerts (in facet).
% local2global gives back the global indices.
function [localContours, local2global, uniqueGlobal] = localIndices(varargin)

numContours = numel(varargin);

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

% I need to create a right-handed coordinate system in which normalVector
% is "up", then project vIn onto this perpspace of normalVector.
function vOut = flatten(vIn, normalVector)

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




