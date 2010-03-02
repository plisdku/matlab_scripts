% createPlaneWaveSource  Initialize a TFSFSource object for a Trogdor simulation
%   createPlaneWaveSource(field, polarization, direction, extent, data) sets up
%   a "soft" total-field scattered-field source.  The given data is written to a
%   file and Trogdor loads it at runtime.  The data should be one value per
%   timestep.
%
%   createPlaneWaveSource(field, polarization, direction, extent, formula) sets
%   up a source using a formula string.
%
%   Further arguments specify omitted sides as axis-aligned unit vectors.
%
%   Example: make a sinusoidal source traveling in the +x direction in a box
%   with Ez polarization, and omit the +x side.
%   
%   s = createPlaneWaveSource('electric', [0 0 1], [0 0 1], [0 0 0 0 0 0], ...
%      'sin(n/20)', [1 0 0]);
%   g = createGrid(..., s);
%
%   See also: createSource, createOneFieldInput, createLink, createGrid
%
%   version 4.5
%   July 29, 2008
function source = createPlaneWaveSource(field, polarization, direction, ...
    dimensions, dataOrFormula, varargin)

source.type = 'TFSFSource';
source.class = 'PlaneWave';
source.fieldName = field;
source.polarization = polarization;
source.direction = direction;
source.dimensions = dimensions;

if ischar(dataOrFormula)
    source.formula = dataOrFormula;
else
    source.data = dataOrFormula;
end

if (nargin >= 6)
    omitSides = varargin;
else
    omitSides = {};
end

source.omitSides = omitSides;


