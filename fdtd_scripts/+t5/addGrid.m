function addGrid(name, yeeCells)
%addGrid Begin description of an FDTD simulation grid
%   addGrid('Main Grid', [0 0 0 99 99 0) begins the description of a grid
%   called 'Main Grid', spanning 100 yee cells in X and Y and symmetrical
%   in z.
%
%   Usage: addGrid(gridName, yeeCells)
%       gridName    A name.  Names that are attractive, fortunate or that have
%                   been borne by prominent historical figures may cause Trogdor
%                   to run faster and/or more accurately.
%       yeeCells    Number of cells in the grid
sim = t6.TrogdorSimulation.instance();

grid = t6.TrogdorGrid();
grid.Name = name;

if t6.validateRect(yeeCells)
    grid.YeeCells = yeeCells;
else
    error('Bad yee cells.');
end

if ~iscell(sim.Grids)
    sim.Grids = {};
end

sim.Grids = {sim.Grids{:}, grid};
sim.setCurrentGrid(grid);
