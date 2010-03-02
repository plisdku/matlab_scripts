% createLink    Initialize a Total-Field/Scattered-Field link for Trogdor
%   createLink(srcGridName, srcExtent, destExtent) will set up a block in the
%   source Grid to provide field values in a region of the destination grid.
%   The srcExtent and destExtent may be the same size, or the srcExtent can
%   have one or two singular dimensions.  The Link object must become a child
%   object of the destination Grid.
%
%   Further arguments specify omitted sides as axis-aligned unit vectors.
%
%   Example: use a 1D source Grid and a 3D destination Grid, and omit the -z
%   side
%
%   srcLink = createLink('Source', [0 0 0 0 0 25], [-25 -25 0 25 25 25], ...
%       [0 0 -1]);
%   mainGrid = createGrid('Main', [10 10 10 10 10 10], ..., srcLink);
%   srcGrid = createGrid('Source', [0 0 10 0 0 10], ...);
%
%   See also: createSource, createPlaneWaveSource, createOneFieldInput, createGrid
%
%   version 4.5
%   July 29, 2008
function link = createLink(srcGridName, srcDimensions, destDimensions, varargin)

if (nargin >= 4)
    omitSides = varargin;
else
    omitSides = {};
end

link.type = 'Link';
link.src = srcGridName;
link.srcDimensions = srcDimensions;
link.destDimensions = destDimensions;
link.omitSides = omitSides;

