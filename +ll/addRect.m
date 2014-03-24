function addRect(varargin)

X.Bounds = [];
X.Material = 1;
X.Conductivity = 0;
X.HMax = '';

X = parseargs(X, varargin{:});

global LL_MODEL;

theMesh = t6.model.Rect(@(p) X.Bounds);

%m = theMesh.meshes(X.Parameters);
%v = m{1}.patchVertices;
%f = m{1}.faces;
%jac = m{1}.jacobian;

meshStruct = struct('material', X.Material, ...
    'bounds', X.Bounds, ...
    'hmax', X.HMax); %, ...
%    'faces', f, 'vertices', v);

LL_MODEL.meshes{numel(LL_MODEL.meshes)+1} = meshStruct;

