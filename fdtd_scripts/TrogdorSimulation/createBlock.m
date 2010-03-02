% createBlock   Initializes a Block object for Trogdor simulation assembly
%   createBlock(dimensions, material) returns a Block object that fills the 
%   given dimensions with one material.
%
%   createBlock(dimensions, material, fillStyle) where fillStyle is either
%   'PECStyle' (the default) or 'PMCStyle' will ensure the appropriate boundary
%   conditions when the Block is written to the simulation grid.
%
%   Example:
%
%   myBlock = createBlock(...);
%   a = createAssembly(..., myBlock);
%
%   See also: createAssembly, createEllipsoid, createKeyImage, createHeightMap
%
%   version 4.5
%   July 29, 2008
function block = createBlock(dimensions, material, varargin)

if (nargin > 2)
    fillStyle = varargin{1};
    assert(strcmp(fillStyle, 'PECStyle') || strcmp(fillStyle, 'PMCStyle'));
else
    fillStyle = 'PECStyle';
end

block.type = 'Block';
block.dimensions = dimensions;
block.material = material;
block.fillStyle = fillStyle;
