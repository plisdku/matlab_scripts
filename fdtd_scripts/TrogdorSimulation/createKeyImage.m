% createKeyImage   Initialize a KeyImage for a Trogdor simulation
%   Usage: createKeyImage(extent, image, rowDirection, colDirection, ...)
%
%   Create a KeyImage oriented such that the rows of the image are oriented
%   along rowDirection and the columns are oriented along colDirection.  The
%   image may be greyscale (NxM) or color (NxMx3) and either floating-point or
%   integer.  Pixels in the image correspond to materials in the simulation.
%   Unused colors in the image will not be filled-in in Trogdor (no change will
%   be made to underlying materials).  Material tags are specified in further
%   arguments as cell arrays, e.g. {'Air', 
%
%   Example: make a funky 2D simulation with grayscale key image
%
%   keyArray = zeros(50, 50);  % zeros will be kryptonite goo
%   keyArray(3:6, 4:7) = 1;  % this part will be PEC
%   keyArray(20:21, 2:40) = 0.5;  % this part will be gold
%   keyImg = createKeyImage([-20 -20 0 29 29 0], keyArray, [1 0 0], [0 -1 0], ...
%       {0, 'Goo', 'PMCStyle'}, {1, 'PEC'}, {0.5, 'Gold'});
%   
%   mainAssembly = createAssembly(..., keyImg);
%
%   See also: createBlock, createEllipsoid, createHeightMap, createAssembly
%
%   version 4.5
%   July 29, 2008
function keyimage = createKeyImage(dimensions, image, imageRow, imageCol, ...
    varargin)

keyimage.type = 'KeyImage';
keyimage.tags = {};
keyimage.materials = {};
keyimage.fillstyle = {};
keyimage.row = imageRow;
keyimage.col = imageCol;
keyimage.dimensions = dimensions;
keyimage.image = image;

for nn = 1:nargin-4
    try
        keyimage.tags = {keyimage.tags{:}, varargin{nn}{1}};
        keyimage.materials = {keyimage.materials{:}, varargin{nn}{2}};
        if length(varargin{nn}) >= 3
            keyimage.fillstyle = {keyimage.fillstyle{:}, varargin{nn}{3}};
        else
            keyimage.fillstyle = {keyimage.fillstyle{:}, 'PECStyle'};
        end
    catch
        warning('Crap.');
    end
end
