function addMesh(mesh)
% Store the vertices and faces of a mesh in the simulation data structure.
% Vertices which touch or reach into the PML will be extended to the outer
% boundary of the PML.  This will distort the geometry of structures which
% do not intersect the PML at right angles.
%
% Example:
%
% newMaterial('Air', 'Numerator', 1, 'Denominator', 1);
% r = t6.model.Rect(@(p) [0 0 0 1 1 1], 'Permittivity', 'Air');
% 
% [finish this later]
% 
% 

import t6.*

sim = simulation();

if iscell(mesh)
    if numel(mesh) > 1
        error('Please add one mesh at a time');
    end
    
    mesh = mesh{1}; % now it's not a cell array!  :-D
end

if isa(mesh, 't6.model.Mesh')
    if ~isempty(sim.CurrentGrid.ParameterizedMeshes)
        error('Mixing meshes and parameterized meshes');
    end
    sim.CurrentGrid.Meshes{end+1} = mesh;
elseif isa(mesh, 't6.model.Node')
    if ~isempty(sim.CurrentGrid.Meshes)
        error('Mixing meshes and parameterized meshes');
    end
    sim.CurrentGrid.ParameterizedMeshes{end+1} = mesh;
end


%{
if nargin == 1
    if iscell(varargin{1})
        if numel(varargin{1}) > 1
            error('Please add one mesh at a time');
        else
            theMesh = varargin{1}{1};
        end
    else
        theMesh = varargin{1};
    end
    X.Vertices = theMesh.patchVertices;
    X.Faces = theMesh.faces;
    X.FreeDirections = theMesh.freeDirections();
    X.Permittivity = theMesh.permittivity;
    X.Permeability = theMesh.permeability;
else
    X.Vertices = [];
    X.Faces = [];
    X.FreeDirections = [];
    X.Permittivity = '';
    X.Permeability = '';
    X = parseargs(X, varargin{:});
end

obj = struct;
obj.type = 'Mesh';

if ~isempty(X.Permittivity)
    obj.permittivity = X.Permittivity;
end

if ~isempty(X.Permeability)
    obj.permeability = X.Permeability;
end

obj.vertices = extendIntoPML(X.Vertices);

obj.faces = X.Faces-1;

if isempty(X.FreeDirections)
    obj.vertexFreeDirections = zeros(size(obj.vertices));
else
    obj.vertexFreeDirections = X.FreeDirections;
end

%patch('Vertices', obj.vertices, 'Faces', obj.faces+1, 'FaceColor', 'g', ...
%    'FaceAlpha', 0.3);
%pause(0.01);

sim.CurrentGrid.Assembly{end+1} = obj;

%}



%{
function v = extendIntoPML(vertices)

sim = t6.simulation();

if size(vertices, 2) ~= 3
    error('Vertex array must be Nx3');
end

v = vertices;

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
%}






