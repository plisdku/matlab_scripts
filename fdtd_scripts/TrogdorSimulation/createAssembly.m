% createAssembly    Create Assembly block of Trogdor parameter file
%   createAssembly(...) takes Blocks, Ellipsoids, HeightMaps and KeyImages from
%   the corresponding create functions and puts them into a list to be included
%   in a Grid object.  Example:
%
%   b = createBlock(...);
%   s = createEllipsoid(...);
%   assembly = createAssembly(b, s);
%   g = createGrid(..., assembly);
%
%   See also: createBlock, createEllipsoid, createHeightMap, createKeyImage
%
%   version 4.5
%   July 29, 2008
function assembly = createAssembly(varargin)

assembly.type = 'Assembly';
assembly.children = varargin;


