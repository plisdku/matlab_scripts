function meshToSTL(mesh, fname)

if ~isa(mesh, 't6.model.Mesh')
    mesh = mesh.meshes;
end

if iscell(mesh)
    mesh = mesh{1};
end

assert(isa(mesh, 't6.model.Mesh'));

faces = mesh.faces;
pv = mesh.patchVertices;

writeSTL(pv, faces, fname);

%{

fh = fopen(fname, 'w');

fprintf(fh, 'solid s\n');

for ff = 1:size(faces, 1)
    
    v0 = pv(faces(ff,1),:);
    v1 = pv(faces(ff,2),:);
    v2 = pv(faces(ff,3),:);
    
    normalVector = cross(v1-v0, v2-v0);
    nv = normalVector/norm(normalVector);
    
    fprintf(fh, 'facet normal %i %i %i\n', nv(1), nv(2), nv(3));
    fprintf(fh, '\touter loop\n');
    
    fprintf(fh, '\t\tvertex %i %i %i\n', v0(1), v0(2), v0(3));
    fprintf(fh, '\t\tvertex %i %i %i\n', v1(1), v1(2), v1(3));
    fprintf(fh, '\t\tvertex %i %i %i\n', v2(1), v2(2), v2(3));
    
    fprintf(fh, '\tendloop\nendfacet\n');
    
end

fprintf(fh, 'endsolid s\n');

%}