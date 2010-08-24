function addLink(sourceGridName, sourceYeeCells, destYeeCells, varargin)
%addLink Add a TFSF boundary to the current grid, sourcing fields from another grid
%   addLink('AuxGrid', [0 0 0 100 0 0], [0 0 0 100 100 100]) will place a TFSF
%       boundary in the current grid that maps the 1D region [0 0 0 100 0 0] of
%       the source grid 'AuxGrid' to the 3D region [0 0 0 100 100 100] of the
%       current grid.  Any mapping of N-D to M-D is permitted where N <= M.
%
%   Usage: addLink(sourceGridName, sourceCells, destCells, named parameters)
%   sourceGridName must be the name of a grid that is added to the simulation
%   (either before or after the current grid) using addGrid.
%
%   Named parameters:
%       OmitSide    An axis-aligned unit vector or cell array of axis-aligned
%                   unit vectors.  The unit vector [1 0 0] specifies the +X
%                   side of the TFSF boundary and instructs Trogdor to not
%                   add or subtract electromagnetic fields on that side.
%                   (default: empty cell array)
%
%   Example:
%
%   addLink('AuxGrid', [0 0 0 100 0 0], [0 0 0 100 100 100], 'OmitSide', ...
%       {[1 0 0], [0 -1 0]});
%
grid = t5.TrogdorSimulation.instance().currentGrid();

if ~t5.validateRect(sourceYeeCells)
    error('Invalid source rectangle.');
end

if ~t5.validateRect(destYeeCells)
    error('Invalid destination rectangle.');
end

obj.type = 'Link';
obj.sourceGrid = sourceGridName;
obj.sourceYeeCells = sourceYeeCells;
obj.destYeeCells = destYeeCells;
obj.omitSides = {};

if length(varargin) == 2
    if ~strcmp(varargin{1}, 'OmitSide')
        error('Only optional attribute is OmitSide');
    end
    
    if ~iscell(varargin{2})
        obj.omitSides = {varargin{2}};
    else
        obj.omitSides = varargin{2};
    end
elseif length(varargin) > 2
    error('Multiple omitted sides must be specified in a cell array.');
end

grid.Links = {grid.Links{:}, obj};

