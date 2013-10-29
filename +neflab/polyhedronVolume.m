function vol = polyhedronVolume(vertices, faces)

% no idea why I bother doing this translation.  :-/  does it get more
% accurate?
%center = mean(vertices);
%displacement = bsxfun(@plus, vertices, -center);

vol = 0;

numFaces = size(faces,1);

for ff = 1:numFaces
    
    faceVerts = vertices(faces(ff,:),:);
    
    vol = vol + det(faceVerts);
    
end

vol = 0.5*vol;
