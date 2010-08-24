function newMaterial(name, zNumerator, zDenominator)
%newMaterial Declare a new permittivity or permeability type
%   newMaterial('Air', 1, 1) creates a material identical to free space
%   newDielectric('Water', 'epsr', 1.7689) declares a dielectric material
%       with epsr = 1.7689 and the default value mur = 1.
%
%   Usage: newMaterial(name, zNumerator, zDenominator)
%       name            a unique identifying string
%       zNumerator      polynomial coefficients in numerator
%       zDenominator    polynomial coefficients in denominator
%

sim = t6.TrogdorSimulation.instance();

%{
X.numerator = [];
X.denominator = [];
X = parseargs(X, varargin{:});

material = struct;
material.name = name;
material.numerator = X.numerator;
material.denominator = X.denominator;
%}
material = struct('name', name, 'numerator', zNumerator, ...
    'denominator', zDenominator);

sim.Materials = {sim.Materials{:}, material};