function newDielectric(name, varargin)
%newDielectric Declare a new lossless, nondispersive dielectric material
%   newDielectric('Air') declares a dielectric material with the default
%       parameters epsr = 1 and mur = 1.
%   newDielectric('Water', 'epsr', 1.7689) declares a dielectric material
%       with epsr = 1.7689 and the default value mur = 1.
%
%   Named parameters:
%       epsr    relative electric permittivity (default 1.0)
%       mur     relative magnetic permeability (default 1.0)
%

sim = t5.TrogdorSimulation.instance();

X.epsr = 1.0;
X.mur = 1.0;
X = parseargs(X, varargin{:});

material = struct;
material.name = name;
material.model = 'StaticDielectric';
material.epsr = X.epsr;
material.mur = X.mur;

sim.Materials = {sim.Materials{:}, material};
