function [outVertices, outFaces, outParents] = trianglePaternity(...
    vChild, fChild, vParents, fParents)
% trianglePaternity     Identify intersections of a triangle with several parent triangles
%
% [v, f, parents] = trianglePaternity(vChild, fChild, vParents, fParents) returns a
% triangulation of the child triangle, labeled by the index of intersecting parent
% triangles.
%
% Parts of the child triangle that intersect no parent triangle are excluded from the
% returned triangulation, so it's up to the caller to make sure the child is really
% covered by the parents.

%% 

%% Triangulate all the child and parent vertices at once
%
% Some of these triangles belong to the child triangle.

numChildVertices = size(vChild, 1);

vertices = [vChild; vParents];
faces = [fChild; fParents + numChildVertices];

edgeConstraints = [faces(:,1:2); faces(:, 2:3); faces(:, [3 1])];

parentTris = triangulation(fParents, vParents);
DT = delaunayTriangulation(vertices, edgeConstraints);
%%
%patch('Vertices', DT.Points, 'Faces', DT.ConnectivityList, 'FaceColor', 'r');

%% Select tris in the child tri

vv = DT.Points;
ff = DT.ConnectivityList;

%% 

centroid = @(face) mean(vv(face,:));

outFaces = [];
outParents = [];
outVertices = DT.Points;

for nn = 1:size(ff, 1)
    triCenter = centroid(ff(nn,:)); % center of current delaunay tri
    
    
    if inpolygon(triCenter(1), triCenter(2), vChild(:,1), vChild(:,2))
        %isPartOfChild(nn) = 1;
        iParent = parentTris.pointLocation(triCenter);
        %fprintf('Parent is %i', parentTri(nn))
        
        if ~isnan(iParent)
            outFaces = [outFaces; ff(nn,:)];
            outParents = [outParents; iParent];
        end
    end
    
    %clf
    %patch('Vertices', DT.Points, 'Faces', DT.ConnectivityList, 'FaceColor', 'r');
    %hold on
    %plot(triCenter(1), triCenter(2), 'o');
    %pause
end

%%
% At this point, parentTri is
% - 0 for delaunay tris that are not in the child tri
% - nonzero (an index) for delaunay tris in the child AND in a parent
% - NaN for delaunay tris not in a parent
%
% so, the paternity test is done and we can return the delaunay tris
% and the parent tri indices.


%%
end

