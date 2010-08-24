function newLossyDielectric(name, varargin)
%newLossyDielectric Declare a new dielectric material with DC conductivity
%   newLossyDielectric('Air') declares a "lossy" dielectric material with the
%       default parameters epsr = 1, mur = 1, and sigma = 0
%   newLossyDielectric('Conductive Air', 'sigma', 100) declares a dielectric
%       material sigma = 100 siemens/m and the default epsr = 1 and mur = 1.
%
%   Named parameters:
%       epsr    relative electric permittivity (default 1.0)
%       mur     relative magnetic permeability (default 1.0)
%       sigma   DC electrical conductivity (default 0.0)

sim = t5.TrogdorSimulation.instance();

X.epsr = 1.0;
X.mur = 1.0;
X.sigma = 0.0;
X = parseargs(X, varargin{:});

material = struct;
material.name = name;
material.model = 'StaticLossyDielectric';
material.epsr = X.epsr;
material.mur = X.mur;
material.sigma = X.sigma;

sim.Materials = {sim.Materials{:}, material};
