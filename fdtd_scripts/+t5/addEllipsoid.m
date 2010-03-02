function addEllipsoid(materialName, yeeCells, varargin)
%addEllipsoid Write materials into an ellipsoidal region of the FDTD grid
%   addEllipsoid('Gold', [10 10 10 20 20 20]) will write a ball of radius 5.
%   addEllipsoid('Air', [10 10 0 20 20 0]) will write a circle of radius 5 in
%       the x-y plane.
%   
%   Usage: addEllipsoid(materialName, yeeCells, named parameters)
%       materialName    Name of a material specified with, e.g., newDielectric()
%       yeeCells        The bounds of the region in which the ellipsoid is given;
%                       [x0 y0 z0 x1 y1 z1] will paint in all cells (x, y, z)
%                       where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.
%
%   Named parameters:
%       FillStyle       'PECStyle' or 'PMCStyle' (default: 'PECStyle')
%
grid = t5.TrogdorSimulation.instance().currentGrid();

if ~t5.validateRect(yeeCells)
    error('Invalid rectangle.');
end

obj = struct;
obj.type = 'Ellipsoid';
obj.materialName = materialName;
obj.yeeCells = yeeCells;

if nargin > 2
    X.FillStyle = {'PECStyle', 'PMCStyle'};
    X = parseargs(X, varargin{:});
    obj.fillStyle = X.FillStyle;
end

grid.Assembly = {grid.Assembly{:}, obj};

