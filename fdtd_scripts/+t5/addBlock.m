function addBlock(materialName, yeeCells, varargin)
%addBlock Fill a rectangular region of the FDTD grid with a material
%   addBlock('Gold', [10 0 0 20 20 0]) will fill a region of cells with 'Gold'.
%
%   Usage: addBlock(materialName, yeeCells, named parameters)
%       materialName    Name of a material specified with, e.g., newDielectric()
%       yeeCells        The bounds of the rect to be filled; [x0 y0 z0 x1 y1 z1]
%                       will paint in all cells (x, y, z) where x0 <= x <= x1,
%                       y0 <= y <= y1, z0 <= z <= z1.
%
%   Named parameters:
%       FillStyle       'PECStyle' or 'PMCStyle' (default: 'PECStyle')
grid = t5.TrogdorSimulation.instance().currentGrid();

if ~t5.validateRect(yeeCells)
    error('Invalid rectangle.');
end

obj = struct;
obj.type = 'Block';
obj.materialName = materialName;
obj.yeeCells = yeeCells;

if nargin > 2
    X.FillStyle = {'PECStyle', 'PMCStyle'};
    X = parseargs(X, varargin{:});
    obj.fillStyle = X.FillStyle;
end

grid.Assembly = {grid.Assembly{:}, obj};

