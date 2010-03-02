% getOutputFrame        Read one frame of data from an open output file
%   [data, count] = getOutputFrame(dataFID, dimensions) reads the next
%   frame from the specified file, assuming that each frame has the given
%   dimensions.  FID is a file identifier from FOPEN, but it may be more
%   convenient to get FID and dimensions both from openOutputFile.  The
%   dimensions of data will be as given (dimensions) and count is the
%   number of float32s read from the file.  At the end of the file, count
%   will be set to zero.
%
%   version 4.5
%   July 29, 2008

function [data, count] = getOutputFrame(dataFID, dimensions)


frameVolume = prod(dimensions);

data = zeros(dimensions);

count = 0;

if (feof(dataFID))
    count = 0;
end

try
    [dataline, cownt] = fread(dataFID, frameVolume, 'float32=>float32');
    if (cownt == frameVolume)
        count = cownt;
        data = reshape(dataline, dimensions);
    else
        count = 0;
    end
catch
    disp(lasterr);
    count = 0;
end
    
end

