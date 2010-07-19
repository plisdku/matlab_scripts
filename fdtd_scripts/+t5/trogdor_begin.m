function trogdor_begin(dxyz, dt, numTimesteps)
%trogdor_begin Begin Trogdor simulation description
%   trogdor_begin([dx dy dz], dt, numTimesteps) declares a simulation where each
%   Yee cell has dimensions [dx dy dz] meters and each timestep is dt seconds.
%   This must be the first line of a Trogdor simulation description.
t6.TrogdorSimulation.clear();
sim = t6.TrogdorSimulation.instance();

sim.Dxyz = dxyz;
sim.Dt = dt;
sim.NumT = numTimesteps;
sim.Grids = [];
sim.CurrentGrid = [];
