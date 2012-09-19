function addGrid(name, yeeCells, varargin)
%addGrid Begin description of an FDTD simulation grid
%   addGrid('Main Grid', [0 0 0 99 99 0]) begins the description of a grid
%   called 'Main Grid', spanning 100 yee cells in X and Y and symmetrical
%   in z.
%
%   Usage: addGrid(gridName, yeeCells)
%       gridName    A name.  Names that are attractive, fortunate or that have
%                   been borne by prominent historical figures may cause Trogdor
%                   to run faster and/or more accurately.
%       yeeCells    Number of cells in the grid

import t6.*;

sim = simulation();

grid = t6.TrogdorGrid();
grid.Name = name;

if validateRect(yeeCells)
    grid.YeeCells = yeeCells;
else
    error('Bad yee cells.');
end

if numel(varargin) > 0
    grid.Origin = varargin{1};
end

if numel(varargin) > 1
    grid.PML = varargin{2};
end

if ~iscell(sim.Grids)
    sim.Grids = {};
end

sim.Grids{end+1} = grid;
sim.setCurrentGrid(grid);


