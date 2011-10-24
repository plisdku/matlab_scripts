function setBackground(varargin)
%setBackground Set a default permittivity and/or permeability for regions
%   of the simulation where the material properties are not otherwise
%   defined
%
%   setBackground('Permittivity', 'Air') will set only the background
%       permittivity explicitly, but the permeability must be separately
%       provided in the entire space
%   setBackground('Permeability', 'Iron') will set only the background
%       permeability
%   setBackground('Permittivity', 'Chocolate', 'Permeability', 'Coffee')
%       will begin preparing a Mocha (run FDTD to complete the drink)
%
%   It is an error to provide more than one background material for a grid.
%   To set both permittivity and permeability, make one call to
%   setBackground().

grid = t6.TrogdorSimulation.instance().currentGrid();

obj = struct;
obj.type = 'Background';

X.Permittivity = '';
X.Permeability = '';
X = parseargs(X, varargin{:});

if ~isempty(X.Permittivity)
    obj.permittivity = X.Permittivity;
end

if ~isempty(X.Permeability)
    obj.permeability = X.Permeability;
end

grid.Assembly = {grid.Assembly{:}, obj};