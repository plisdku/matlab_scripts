function addMesh(varargin)
% addMesh    

X.Mesh = [];
X.Material = 1;
X.Voltage = [];
X.SurfaceCharge = [];
X.Conductivity = 0;
X.Parameters = [];
X.HMax = '';
X.HMin = '';
X.HGrad = '';
X.Exclude = 0;

X = parseargs(X, varargin{:});

global LL_MODEL;

m = X.Mesh.meshes(X.Parameters);
v = m{1}.patchVertices;
f = m{1}.faces;
jac = m{1}.jacobian;

if ~isempty(X.Voltage) && ~isempty(X.SurfaceCharge)
    error(['Both Voltage and SurfaceCharge specified (this is ' ...
        'Dirichlet and Neumann)']);
end

%numMeshes = numel(LL_MODEL.meshes);
%fname = sprintf('mesh_%i.stl', numMeshes+1);
%ll.meshToSTL(X.Mesh, fname);

meshStruct = struct('material', X.Material, ...
    'voltage', X.Voltage, 'surfacecharge', X.SurfaceCharge, ...
    'vertices', v, 'faces', f, 'jacobian', jac, ...
    'hmax', X.HMax, 'hgrad', X.HGrad, 'hmin', X.HMin, 'exclude', X.Exclude);

LL_MODEL.meshes{numel(LL_MODEL.meshes)+1} = meshStruct;

