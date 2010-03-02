function addGrid(name)
%addGrid Begin description of an FDTD simulation grid
%   addGrid('Main Grid') begins the description of a grid called 'Main Grid'.
%
%   Usage: addGrid(gridName)
%       gridName    A name.  Names that are attractive, fortunate or that have
%                   been borne by prominent historical figures may cause Trogdor
%                   to run faster and/or more accurately.
sim = t5.TrogdorSimulation.instance();

grid = t5.TrogdorGrid();
grid.Name = name;

if ~iscell(sim.Grids)
    sim.Grids = {};
end

sim.Grids = {sim.Grids{:}, grid};
sim.setCurrentGrid(grid);
