% getOutputFrames       Extract time sequence of output data from FDTD sim.
%   [data, count] = getOutputFrames(FID, dimensions) reads all the frames
%   from the specified file, assuming that each frame has the given
%   dimensions.  FID is a file identifier from FOPEN, but it may be more
%   convenient to get FID and dimensions both from openOutputFile.  The
%   dimensions of data will be [dimensions, numFrames] and count is the
%   number of float32s read from the file.
%
%   version 4.5
%   July 29, 2008

function [data, numFrames] = getOutputFrames(dataFID, dimensions)



[data, count] = fread(dataFID, inf, 'float32=>float32');

%   We must trim off any extraneous incomplete data.  Calculate 
%   how many complete frames (timesteps) are contained in the data
%   and truncate the data before reshaping.
frameSize = prod(dimensions);
numFrames = floor(count / frameSize);

if (numFrames*frameSize ~= length(data(:)))
    data = data(1:(numFrames*frameSize));
end
data = reshape(data, [dimensions, numFrames]);


