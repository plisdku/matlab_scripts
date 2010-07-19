function newDrude(name, varargin)
%newDrude Declare a new Drude metal
%   newDrude('Gold', 'epsinf', 12.99, 'omegap', 4e15, 'tauc', 9e-15) declares
%       a Drude-model material that approximates the behavior of gold in the
%       visible domain (you are advised to carefully choose the parameters
%       yourself!).
%
%   Named parameters:
%       epsinf  infinite-frequency relative electric permittivity (default 1.0)
%       mur     relative magnetic permeability (default 1.0)
%       omegap  plasmon frequency (default 0.0)
%       tauc    relaxation time (default 0.0)
sim = t6.TrogdorSimulation.instance();

X.epsinf = 1;
X.omegap = 0;
X.tauc = 0;
X = parseargs(X, varargin{:});

material = struct;
material.name = name;
material.model = 'DrudeMetal1';
material.epsinf = X.epsinf;
material.omegap = X.omegap;
material.tauc = X.tauc;

sim.Materials = {sim.Materials{:}, material};
