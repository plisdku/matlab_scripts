function writePolyhedron(vertices, faces, fid)

if nargin < 3
    fid = 1;
end

numVertices = size(vertices,1);
numFaces = size(faces,1);

fprintf(fid, '%i %i\n', numVertices, numFaces);

for vv = 1:size(vertices,1)
    fprintf(fid, '%2.8f %2.8f %2.8f\n', ...
        vertices(vv,1), vertices(vv,2), vertices(vv,3));
end

for ff = 1:size(faces,1)
    
    faceForC = faces(ff,:) - 1;
    
    fprintf(fid, '%i %i %i\n', faceForC(1), faceForC(2), faceForC(3));
end