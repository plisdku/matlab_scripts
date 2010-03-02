% createHeightMap   Initialize a HeightMap for a Trogdor simulation
%   Usage: createHeightMap(extent, image, material, rowDirection, colDirection,
%       upDirection, ...)
%
%   Create a heightmap oriented such that the rows of the image are oriented
%   along rowDirection, the columns are oriented along colDirection and the
%   height values increase along upDirection.  It is optional to provide a
%   fillStyle argument.
%
%   The image parameter should take on values from zero to 1 (if floating-point)
%   or zero to 255 (if integer).  The extent of the image should be one cell
%   smaller than size(image), otherwise Trogdor's behavior is undefined.
%
%   EXAMPLE: make a heightmap using peaks() ("Mount Matlab")
%
%   peaksheights = peaks(61);
%   heightmap = peaksheights - min(peaksheights(:));  % lowest value is zero
%   heightmap = heightmap ./ max(heightmap(:)); % scale between zero and 1
%
%   PECheights = createHeightMap([-30 -30 -30 30 30 0], heightmap, 'PEC', ...
%    [1 0 0], [0 1 0], [0 0 1]);
%
%   PMCheights = createHeightMap([-30 -30 0 30 30 30], heightmap, 'PEC', ...
%    [1 0 0], [0 1 0], [0 0 -1], 'PMCStyle');
%    
%   mainAssembly = createAssembly(PECheights, PMCheights);
%   
%   See also: createBlock, createEllipsoid, createKeyImage, createAssembly
%
%   version 4.5
%   July 29, 2008
function heightmap = createHeightMap(dimensions, image, material, ...
    imageRowAxis, imageColAxis, imageUpAxis, varargin)

if (nargin >= 7)
    fillStyle = varargin{1};
else
    fillStyle = 'PECStyle';
end

heightmap.type = 'HeightMap';
heightmap.dimensions = dimensions;
heightmap.image = image;
heightmap.material = material;
heightmap.row = imageRowAxis;
heightmap.col = imageColAxis;
heightmap.up = imageUpAxis;
heightmap.fillStyle = fillStyle;


