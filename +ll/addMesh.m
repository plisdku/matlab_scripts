function addMesh(varargin)

X.Mesh = [];
X.Material = 1;
X.Conductivity = 0;
X.Parameters = [];
X.HMax = '';

X = parseargs(X, varargin{:});

global LL_MODEL;

m = X.Mesh.meshes(X.Parameters);
v = m{1}.patchVertices;
f = m{1}.faces;
jac = m{1}.jacobian;

%numMeshes = numel(LL_MODEL.meshes);
%fname = sprintf('mesh_%i.stl', numMeshes+1);
%ll.meshToSTL(X.Mesh, fname);

meshStruct = struct('material', X.Material, ...
    'vertices', v, 'faces', f, 'jacobian', jac, ...
    'hmax', X.HMax);

LL_MODEL.meshes{numel(LL_MODEL.meshes)+1} = meshStruct;

