function addSTLMesh(varargin)

X.Mesh = [];
X.Material = 1;
X.Conductivity = 0;

X = parseargs(X, varargin{:});

global LL_MODEL;

numMeshes = numel(LL_MODEL.meshes);

fname = sprintf('mesh_%i.stl', numMeshes+1);

ll.meshToSTL(X.Mesh, fname);

meshStruct = struct('fname', fname, ...
    'material', X.Material);

LL_MODEL.meshes{numel(LL_MODEL.meshes)+1} = meshStruct;