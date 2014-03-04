function material(varargin)

X.Permittivity = 1;
X.Permeability = 1;
X.Conductivity = 0;
X.Name = '';

X = parseargs(X, varargin{:});

newMat = struct('name', X.Name, ...
    'epsr', X.Permittivity, ...
    'mur', X.Permeability, ...
    'sigma', X.Conductivity);

global LL_MODEL;

LL_MODEL.materials{numel(LL_MODEL.materials)+1} = newMat;
