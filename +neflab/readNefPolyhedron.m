function [vertices, faces] = readNefPolyhedron(fid)

l = fgetl(fid);
[numVertices count] = sscanf(l, 'numVertices %i');
assert(count == 1);

vertices = zeros(numVertices,3);

%fprintf('%i vertices...\n', numVertices);
for vv = 1:numVertices
    l = fgetl(fid);
    [a, count] = sscanf(l, '%f %f %f %f');
    assert(count == 4);
    
    vertices(vv,:) = a(2:4);
end

%disp(vertices);

l = fgetl(fid);
[numFacets count] = sscanf(l, 'numFacets %i');
assert(count == 1);

%fprintf('%i facets...\n', numFacets);

faces = [];
for ff = 1:numFacets
    
    l = fgetl(fid);
    [numContours count] = sscanf(l, '\tnumContours %i');
    assert(count == 1);
    
    %fprintf('\tFacet %i has %i contours\n', ff, numContours);
    
    contours = cell(numContours,1);
    for cc = 1:numContours
        l = fgetl(fid);
        [numContourVerts count] = sscanf(l, '\tnumVertices %i');
        assert(count == 1);
        
        %fprintf('\t\tContour %i has %i vertices\n', cc, numContourVerts);
        
        loopVerts = zeros(numContourVerts,1);
        l = fgetl(fid);
        
        loopVerts = sscanf(l, '%i');
        assert(numel(loopVerts) == numContourVerts);
        
        loopVerts = loopVerts + 1; % bump from C++ to Matlab indexing
        %disp(loopVerts')
        
        contours{cc} = loopVerts;
        
        %disp(vertices)
        %figure(1); clf
        %patch('Vertices', vertices, 'Faces', loopVerts', 'FaceColor', 'r');
        %view(3)
        %axis vis3d
    end
    
    % To perform the triangulation I need to translate all the vertex
    % indices in the contour array.  Hm.
    
    tris = neflab.triangulate(vertices, contours{:});
    faces = [faces; tris];
    %disp(tris)
    %pause
end








































