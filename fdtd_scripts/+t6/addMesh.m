function addMesh(varargin)
% Store the vertices and faces of a mesh in the simulation data structure.
% Vertices which touch or reach into the PML will be extended to the outer
% boundary of the PML.  This will distort the geometry of structures which
% do not intersect the PML at right angles.

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

obj.vertices = extendIntoPML(X.Vertices);

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



function v = extendIntoPML(vertices)

if size(vertices, 2) ~= 3
    error('Vertex array must be Nx3');
end

v = vertices;

sim = t6.simulation();
outerBounds = sim.OuterBounds;
innerBounds = sim.NonPMLBounds;

for xyz = 1:3
if innerBounds(xyz) ~= innerBounds(xyz+3)
    iFloor = vertices(:,xyz) <= innerBounds(xyz);
    iCeil = vertices(:,xyz) >= innerBounds(xyz+3);
    
    v(iFloor,xyz) = outerBounds(xyz);
    v(iCeil,xyz) = outerBounds(xyz+3);
end
end







