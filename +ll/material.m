function material(varargin)

X.Permittivity = 1;
X.Permeability = 1;
X.Conductivity = 0;
X.Name = '';

X = parseargs(X, varargin{:});

if imag(X.Permittivity) > 0
    warning('Permittivity has positive imaginary part: gain in COMSOL!!');
end

newMat = struct('name', X.Name, ...
    'epsr', X.Permittivity, ...
    'mur', X.Permeability, ...
    'sigma', X.Conductivity);

global LL_MODEL;

LL_MODEL.materials{numel(LL_MODEL.materials)+1} = newMat;
