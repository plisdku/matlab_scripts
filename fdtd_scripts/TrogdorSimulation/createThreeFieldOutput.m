% createThreeFieldOutput     Initialize an output for a Trogdor simulation
%   createThreeFieldOutput(fieldName, dimensions, filename) sets up an output
%   of the given field ('electric' or 'magnetic') in the region specified by
%   dimensions, to the given filename.
%
%   createThreeFieldOutput(fieldName, dimensions, filename, period) will save
%   data at timesteps [period, 2*period, 3*period, ...].
%
%   Example:
%
%   myOutput = createThreeFieldOutput('magnetic', [x1 y1 z1 x2 y2 z2], 'foo');
%   g = createGrid(..., myOutput);
%
%   ThreeFieldOutput saves the three electric or magnetic field components
%   exactly as they are in the simulation without resampling to account for
%   the Yee cell structure.
%
%   See also: createOneFieldOutput, createColocatedOutput, createOneFieldInput
%
%   version 4.5
%   July 29, 2008
function output = createThreeFieldOutput(fieldName, dimensions, filename, varargin)

if (nargin >= 4)
    output.period = varargin{1};
else
    output.period = 1;
end

if (nargin >= 5)
    output.stride = varargin{2};
else
    output.stride = [1 1 1];
end

output.type = 'Output';
output.class = 'ThreeFieldOutput';
output.field = fieldName;
output.dimensions = dimensions;
output.filename = filename;
