% createMaterial    Initialize a Material for a Trogdor simulation
%   createMaterial(name, className) creates a Material object with no material
%   parameters.
%
%   createMaterial(name, className, parameter1, ...) creates a Material object
%   with material parameters.  Each parameter is a cell array pair of strings.
%   Note that numerical parameters must be passed as strings here.
%
%   Example:
%
%   drudeMaterial = createMaterial('Gold', 'DrudeMetalModel', {'epsinf', ...
%       '12.9898'}, {'omegap', '4e15'}, {'tauc', '9e-15'});
%
%   See also: createDrudeMaterial, createPEC, createStaticDielectric, createStaticLossyDielectric
%
%   version 4.5
%   July 29, 2008
function material = createMaterial(name, class, varargin)

material.type = 'Material';
material.name = name;
material.class = class;
material.parameters = varargin;
