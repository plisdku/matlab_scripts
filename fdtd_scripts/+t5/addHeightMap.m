function addHeightMap(materialName, yeeCells, image, imageRow, imageCol, ...
    imageUp, varargin)
%addHeightMap Write materials into the FDTD grid using a heightmap
%   addHeightMap('Gold', [0 0 0 99 99 99], heightMap, [1 0 0], [0 1 0], ...
%       [0 0 1]) will fill all cells (x, y, z) in the region [0 0 0 99 99 99]
%       where 100*heightMap(x, y) <= z with 'Gold'.
%   
%   Usage: addHeightMap(material, yeeCells, heightMap, xDir, yDir, zDir, ...)
%       material    Name of a material specified with, e.g., newDielectric()
%       yeeCells    The bounds of the region in which the heightmap is given;
%                   [x0 y0 z0 x1 y1 z1] will paint in all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.
%       heightMap   An array of size [nx ny nz] taking on values from 0 to 1.
%                   Here nx is the size of the Yee cell region along xDir,
%                   ny is the size of the Yee cell region along yDir, and
%                   nz is the height of the Yee cell region along zDir.  The
%                   lowest value 0 will not fill any cells, and the highest
%                   value 1 will fill all cells.  THIS ARRAY IS IN CARTESIAN
%                   COORDINATES; in particular columns of this array point
%                   along +y by default, as opposed to -y as conventional for
%                   images.  This is to remain consistent with Trogdor's input
%                   and output coordinate system.
%       xDir, yDir, zDir    Axis-aligned unit vectors that determine how the
%                   heightMap is oriented in the FDTD space; xDir and yDir
%                   give transverse directions and zDir is "up".
%
%   Named parameters:
%       FillStyle   either 'PECStyle' or 'PMCStyle' (default: 'PECStyle')
%
grid = t5.TrogdorSimulation.instance().currentGrid();

if ~t5.validateRect(yeeCells)
    error('Invalid rectangle.');
end

if ndims(image) ~= 2
    error('Height map image must be 2D.');
end

if max(image(:)) > 1 || min(image(:)) < 0
    error('Height map values must range from 0 to 1.');
end

obj = struct;
obj.type = 'HeightMap';
obj.materialName = materialName;
obj.yeeCells = yeeCells;
obj.row = t5.makeAxisString(imageRow);
obj.column = t5.makeAxisString(imageCol);
obj.up = t5.makeAxisString(imageUp);
obj.image = image;

if length(varargin) == 2
    if ~strcmp(varargin{1}, 'FillStyle')
        error('Only optional HeightMap attribute is FillStyle.');
    end
    
    if ~strcmp(varargin{2}, 'PECStyle') && ~strcmp(varargin{2}, 'PMCStyle')
        error('Only valid FillStyles are PECStyle and PMCStyle.');
    end
    
    obj.fillStyle = varargin{2};
end

grid.Assembly = {grid.Assembly{:}, obj};

