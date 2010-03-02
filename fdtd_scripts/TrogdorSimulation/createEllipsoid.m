% createEllipsoid   Initializes an Ellipsoid object for Trogdor simulation
%   createEllipsoid(dimensions, material) returns an Ellipsoid object that fills
%   the given dimensions with one material.
%
%   createEllipsoid(dimensions, material, fillStyle) where fillStyle is either
%   'PECStyle' (the default) or 'PMCStyle' will ensure the appropriate boundary
%   conditions when the Block is written to the simulation grid.
%
%   Example:
%
%   mySphere = createEllipsoid([-5 -5 -5 5 5 5], 'Gold');
%   a = createAssembly(mySphere);
%
%   See also: createBlock, createKeyImage, createHeightMap, createAssembly
%
%   version 4.5
%   July 29, 2008
function ball = createEllipsoid(dimensions, material, varargin)

if (nargin == 4)
    fillStyle = varargin{1};
    assert(strcmp(fillStyle, 'PECStyle') || strcmp(fillStyle, 'PMCStyle'));
else
    fillStyle = 'PECStyle';
end

ball.type = 'Ellipsoid';
ball.dimensions = dimensions;
ball.material = material;
ball.fillStyle = fillStyle;

