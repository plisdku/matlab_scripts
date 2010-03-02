% createGrid    Initialize a simulation Grid for a Trogdor simulation
%   createGrid(gridName, PML, ...) collects assembly, input/output/source and 
%   link objects to run an FDTD simulation of a single structure.  Several Grids
%   may be run in one Simulation, which is only useful for setting up TFSF
%   soft sources.  Arguments should include, in arbitrary order, Assembly,
%   Link, Source, Input, and Output objects.
%
%   Example:
%
%   a = createAssembly(...);
%   s = createSource(...);
%   o = createColocatedOutput(...);
%   mainGrid = createGrid('Main Grid', [10 10 0 10 10 0], a, s, o);
%   s = createSimulation([dx dy dz dt], numT, mainGrid);
%
%   See also: createAssembly, createSimulation
%
%   version 4.5
%   July 29, 2008
function grid = createGrid(gridName, PML, varargin)

grid.type ='Grid';
grid.name = gridName;
grid.PML = PML;
grid.pieces = varargin;

