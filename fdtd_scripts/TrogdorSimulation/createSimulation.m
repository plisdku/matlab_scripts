% createSimulation  Initialize a Trogdor simulation
%   createSimulation(dxdydzdt, numTimesteps, grid1, grid2, ...) combines all the
%   necessary child objects for a Trogdor simulation.  This object is ready
%   to be passed to writeSimulation.
%
%   Example: set up a 1000-timestep simulation with a main grid and source grid
%
%   dx = 10e-9;
%   dt = 0.99 * dx / (sqrt(3.0) * 2.99e8);  % 3D grid Courant limit
%   
%   g = createGrid('MainGrid', ...);
%   s = createGrid('SourceGrid', ...);
%   simObject = createSimulation([dx, dx, dx, dt], 1000, g, s);
%   writeSimulation(simObject);
%   
%   See also: createGrid, writeSimulation
%
%   version 4.5
%   July 29, 2008
function simObject = createSimulation(dxdydzdt, numTimesteps, varargin)

dx = dxdydzdt(1);
dy = dxdydzdt(2);
dz = dxdydzdt(3);
dt = dxdydzdt(4);

simObject.type = 'Simulation';
simObject.dx = dx;
simObject.dy = dy;
simObject.dz = dz;
simObject.dt = dt;
simObject.numT = numTimesteps;

simObject.materials = {};
simObject.grids = {};

for (nn = 1:nargin-2)
    switch (varargin{nn}.type)
        case 'Material'
            simObject.materials = {simObject.materials{:}, varargin{nn}};
        case 'Grid'
            simObject.grids = {simObject.grids{:}, varargin{nn}};
        otherwise
            warning(sprintf('Ignoring child with type %s', varargin{nn}.type));
    end
end

