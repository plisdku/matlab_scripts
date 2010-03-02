function addKeyImage(yeeCells, image, imageRow, imageCol, varargin)
%addKeyImage Write materials into the FDTD grid using a color-coded key image.
%   addKeyImage([10 10 0 60 60 20], keyImage, [1 0 0], [0 1 0], ...
%       {[1 0 0], 'Gold'}, {[0 0 0], 'Vacuum'}) will interpret the values of
%       keyImage as materials to be written in the x-y plane and extruded
%       along z, where red (color [1 0 0]) stands for 'Gold' and black ([0 0 0])
%       stands for 'Vacuum'.
%
%   Usage: addKeyImage(yeeCells, keyImage, xDir, yDir, named parameters)
%       yeeCells    a single rectangle [x0 y0 z0 x1 y1 z1]
%       keyImage    a greyscale or RGB image, of size [nx ny] or [nx ny 3]
%                   respectively.  THE PIXELS OF THIS IMAGE ARE INTERPRETED IN
%                   CARTESIAN COORDINATES consistent with Trogdor input and
%                   output data.  Thus keyImage(2, 10) refers to a material at
%                   x = 2 and y = 10.
%       xDir, yDir  Specify the plane of the image in the Trogdor simulation.
%                   Materials will extrude along the third direction.  These
%                   arguments must be perpendicular axis-aligned unit vectors,
%                   e.g. [1 0 0] and [0 -1 0] which interprets the columns of
%                   keyImage as pointing down along y as for an image.
%
%   Extra parameters: the colors of keyImage are interpreted according to tags
%   which are specified as the last arguments to addKeyImage.  Each tag is
%   of the form {color, materialName} or {color, materialName, fillStyle} where:
%       color           scalar or RGB pixel value in keyImage
%       materialName    name of a simulated material, as specified in (e.g.)
%                       newDielectric, newDrude, etc.
%       fillStyle       'PECStyle' (default if unspecified) or 'PMCStyle'
%                       Do not specify fillStyle if you do not know what it
%                       does.
%
%   Pixels in greyscale keyImage data must vary from 0 to 1.
%   Pixels in RGB keyImage data must vary from 0 to 1 in each channel (R, G, B).

grid = t5.TrogdorSimulation.instance().currentGrid();

if ~t5.validateRect(yeeCells)
    error('Invalid rectangle.');
end

if ndims(image) ~= 2 && ndims(image) ~= 3
    error('Key image must be 2D with optional third RGB dimension.');
end

obj = struct;
obj.type = 'KeyImage';
obj.yeeCells = yeeCells;
obj.row = t5.makeAxisString(imageRow);
obj.column = t5.makeAxisString(imageCol);
obj.image = image;
obj.tags = {};

for tt = 1:length(varargin)
    tag = varargin{tt};
    
    if length(tag) < 2
        error('KeyImage tag must be a cell array like {0.5, materialName}.');
    end
    if ~isnumeric(tag{1})
        error('KeyImage tag must have the form {pixel, materialName}.');
    end
    
    if length(tag{1}) == 1
        if tag{1} < 0 || tag{1} > 1
            error('Greyscale tag must be between 0 and 1.');
        end
    elseif length(tag{1}) == 3
        if any(tag{1} < 0) || any(tag{1} > 1)
            error('Color tag RGB values must be between 0 and 1.');
        end
    else
        error('Tag color must be greyscale or RGB.');
    end
    
    if length(tag) == 2
        tagStruct = struct('pixel', tag{1}, 'material', tag{2}, ...
            'fillStyle', 'PECStyle');
    else
        if ~strcmp(tag{3}, 'PECStyle') && ~strcmp(tag{3}, 'PMCStyle')
            error('KeyImage tag fill style must be PECStyle or PMCStyle.');
        end
        tagStruct = struct('pixel', tag{1}, 'material', tag{2}, ...
            'fillStyle', tag{3});
    end
    
    obj.tags = {obj.tags{:}, tagStruct};
end

grid.Assembly = {grid.Assembly{:}, obj};

