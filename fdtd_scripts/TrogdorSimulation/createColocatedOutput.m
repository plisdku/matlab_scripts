% createColocatedOutput     Initialize an output for a Trogdor simulation
%   createColocatedOutput(fieldName, dimensions, filename) sets up an output
%   of the given field ('electric' or 'magnetic') in the region specified by
%   dimensions, to the given filename.
%
%   ColocatedOutput resamples electric and magnetic field components so they are
%   effectively measured at the same location in the Yee cell, i.e. it undoes
%   the interleaving in space.  (The samples are still interleaved in time.)
%   
%   createColocatedOutput(fieldName, dimensions, filename, period) will save
%   data at timesteps [period, 2*period, 3*period, ...].
%
%   createColocatedOutput(fieldName, dimensions, filename, period, stride) where
%   stride is a three-element vector will decimate the saved data in space.
%   This is useful for coarsely sampling a large output region.  stride should
%   be a row vector of integers, e.g. [10 10 1] to coarsely sample in X and Y
%   and finely sample in Z.  Sampling will include the bottom corner of the
%   dimensions, i.e. dimensions(1:3), and will include the upper corner only if
%   the distance from the bottom corner is a multiple of the stride in all
%   dimensions.
%
%   Example:
%
%   myOutput = createColocatedOutput('electric', [x1 y1 z1 x2 y2 z2], 'foo');
%   g = createGrid(..., myOutput);
%
%   See also: createThreeFieldOutput, createOneFieldOutput, createGrid
%
%   version 4.5
%   July 29, 2008
function output = createColocatedOutput(fieldName, dimensions, filename, varargin)

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
output.class = 'ColocatedOutput';
output.field = fieldName;
output.dimensions = dimensions;
output.filename = filename;
