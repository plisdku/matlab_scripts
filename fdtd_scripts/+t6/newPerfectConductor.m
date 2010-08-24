function newPerfectConductor(name)
%newPerfectConductor Declare a new perfectly conducting material
%   newDielectric('Conductor') declares a material suitable for use as PEC or
%       PMC.
%
%   Named parameters:
%       (none)
sim = t6.TrogdorSimulation.instance();

material = struct;
material.name = name;
material.model = 'PerfectConductor';

sim.Materials = {sim.Materials{:}, material};