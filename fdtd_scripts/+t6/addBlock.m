function addBlock(varargin)
%addBlock Fill a rectangular region of the FDTD grid with a material
%   addBlock('Permittivity', 'Gold', 'YeeBounds', [10 0 0 20 20 0]) will
%   fill a region of cells with 'Gold'.
%
%   Usage: addBlock(materialName, yeeBounds, named parameters)
%       materialName    Name of a material specified with, e.g., newDielectric()
%       yeeBounds        The bounds of the rect to be filled; [x0 y0 z0 x1 y1 z1]
%                       will paint in all cells (x, y, z) where x0 <= x <= x1,
%                       y0 <= y <= y1, z0 <= z <= z1.
%
%   Named parameters:
%       FillStyle       'PECStyle' or 'PMCStyle' (default: 'PECStyle')
grid = t6.TrogdorSimulation.instance().currentGrid();

obj = struct;
obj.type = 'Block';

X.YeeBounds = [0 0 0 0 0 0];
X.Permittivity = '';
X.Permeability = '';
X = parseargs(X, varargin{:});

if ~isempty(X.Permittivity)
    obj.permittivity = X.Permittivity;
end

if ~isempty(X.Permeability)
    obj.permeability = X.Permeability;
end

obj.yeeBounds = X.YeeBounds;

if ~t6.validateRect(obj.yeeBounds)
    error('Invalid rectangle.');
end

grid.Assembly = {grid.Assembly{:}, obj};

