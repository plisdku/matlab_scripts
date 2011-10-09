function addMesh(varargin)

grid = t6.TrogdorSimulation.instance().currentGrid();

obj = struct;
obj.type = 'Mesh';

X.Vertices = [];
X.Faces = [];
X.FreeDirections = [];
X.Permittivity = '';
X.Permeability = '';
X = parseargs(X, varargin{:});

if ~isempty(X.Permittivity)
    obj.permittivity = X.Permittivity;
end

if ~isempty(X.Permeability)
    obj.permeability = X.Permeability;
end

obj.vertices = X.Vertices;
obj.faces = X.Faces-1;

if isempty(X.FreeDirections)
    obj.vertexFreeDirections = zeros(size(obj.vertices));
else
    obj.vertexFreeDirections = X.FreeDirections;
end

grid.Assembly = {grid.Assembly{:}, obj};
