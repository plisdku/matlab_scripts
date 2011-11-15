function addMesh(varargin)

grid = t6.TrogdorSimulation.instance().currentGrid();

obj = struct;
obj.type = 'Mesh';

if nargin == 1
    X.Vertices = varargin{1}.patchVertices;
    X.Faces = varargin{1}.faces;
    X.FreeDirections = varargin{1}.freeDirections();
    X.Permittivity = varargin{1}.permittivity;
    X.Permeability = varargin{1}.permeability;
else
    X.Vertices = [];
    X.Faces = [];
    X.FreeDirections = [];
    X.Permittivity = '';
    X.Permeability = '';
    X = parseargs(X, varargin{:});
end

if ~isempty(X.Permittivity)
    obj.permittivity = X.Permittivity;
end

if ~isempty(X.Permeability)
    obj.permeability = X.Permeability;
end

obj.vertices = X.Vertices;

%obj.vertices = bsxfun(@minus, X.Vertices, grid.Origin);
%obj.vertices(1,:) = obj.vertices(1,:);
%obj.vertices(2,:) = obj.vertices(2,:);
%obj.vertices(3,:) = obj.vertices(3,:);
obj.faces = X.Faces-1;

if isempty(X.FreeDirections)
    obj.vertexFreeDirections = zeros(size(obj.vertices));
else
    obj.vertexFreeDirections = X.FreeDirections;
end

grid.Assembly = {grid.Assembly{:}, obj};
