% createOneFieldInput   Initialize OneFieldInput for Trogdor simulation
%   createOneFieldInput(field, extent, filename) returns an Input object to be
%   parented by a Grid object.  Valid field arguments are 'ex', 'ey', 'ez',
%   'hx', 'hy', 'hz'.  The filename should be the prefix (extensionless) name
%   of a valid data (.dat) and specification (.txt) file pair for Trogdor, such
%   as created by a previous Trogdor output.
%
%   Example:
%
%   inp = createOneFieldInput('ex', [10 10 0 20 20 0], 'sineWave');
%
%   g = createGrid('Main Grid', [0 0 0 30 30 0], inp);
%
%   See also: createSource, createPlaneWaveSource, createLink
%
%   version 4.5
%   July 29, 2008
function input = createOneFieldInput(field, dimensions, filename)

input.type = 'OneFieldInput';
input.dimensions = dimensions;
input.filename = filename;
input.class = 'OneFieldInput';
