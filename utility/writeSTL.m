function writeSTL(pv, faces, fname, solidName)
% writeSTL(vertices_Nx3, faces_Mx3, fname)

if nargin < 4
    solidName = 's';
end

% Vectorize some of this.  For big meshes this is annoyingly slow...
iFaces = 1:size(faces,1);
v0s = pv(faces(iFaces,1),:);
v1s = pv(faces(iFaces,2),:);
v2s = pv(faces(iFaces,3),:);

normalVectors = cross(v1s - v0s, v2s - v0s);
nvs = bsxfun(@times, normalVectors, ...
    1./sqrt(sum(normalVectors.*normalVectors, 2)));

fh = fopen(fname, 'w');

fprintf(fh, 'solid %s\n', solidName);

fmt = ['facet normal %i %i %i\n' ...
    '\touter loop\n' ...
    '\t\tvertex %1.15g %1.15g %1.15g\n' ...
    '\t\tvertex %1.15g %1.15g %1.15g\n' ...
    '\t\tvertex %1.15g %1.15g %1.15g\n' ...
    '\tendloop\n' ...
    'endfacet\n'];

allDat = [nvs(:,1), nvs(:,2), nvs(:,3), ...
    v0s(:,1), v0s(:,2), v0s(:,3), ...
    v1s(:,1), v1s(:,2), v1s(:,3), ...
    v2s(:,1), v2s(:,2), v2s(:,3)];

fprintf(fh, fmt, transpose(allDat));

fprintf(fh, 'endsolid %s\n', solidName);

%{

for ff = 1:size(faces, 1)
    
    
    fprintf(fh, 'facet normal %i %i %i\n', nvs(ff,1), nvs(ff,2), nvs(ff,3));
    
    fprintf(fh, '\touter loop\n');
    
    
    fprintf(fh, '\t\tvertex %1.15g %1.15g %1.15g\n', v0s(ff,1), v0(2), v0(3));
    fprintf(fh, '\t\tvertex %1.15g %1.15g %1.15g\n', v1(1), v1(2), v1(3));
    fprintf(fh, '\t\tvertex %1.15g %1.15g %1.15g\n', v2(1), v2(2), v2(3));
    
    fprintf(fh, '\tendloop\nendfacet\n');
    
    %{
    v0 = pv(faces(ff,1),:);
    v1 = pv(faces(ff,2),:);
    v2 = pv(faces(ff,3),:);
    
    normalVector = cross(v1-v0, v2-v0);
    nv = normalVector/norm(normalVector);
    
    fprintf(fh, 'facet normal %i %i %i\n', nv(1), nv(2), nv(3));
    fprintf(fh, '\touter loop\n');
    
    fprintf(fh, '\t\tvertex %1.15g %1.15g %1.15g\n', v0(1), v0(2), v0(3));
    fprintf(fh, '\t\tvertex %1.15g %1.15g %1.15g\n', v1(1), v1(2), v1(3));
    fprintf(fh, '\t\tvertex %1.15g %1.15g %1.15g\n', v2(1), v2(2), v2(3));
    
    fprintf(fh, '\tendloop\nendfacet\n');
    %}
    
end
%}

fclose(fh);


