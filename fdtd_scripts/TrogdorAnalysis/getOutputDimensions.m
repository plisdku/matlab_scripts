% getOutputDimensions   Get size of each frame matrix in given output file
%   dimensions = getOutputDimensions(filePrefix) returns the dimensions of
%   each frame of an output file.  An output that stores a 100 x 100 array
%   of field values each timestep will have dimensions [100, 100].
%
%   See also: getOutputBounds, getOutputPeriod, readOutputFile
%
%   version 4.5
%   July 29, 2008

function dimensions = getOutputDimensions(fileprefix)


dimensions = [-1, -1, -1];

datfile = [fileprefix, '.dat'];
specfile = [fileprefix, '.txt'];


specFID = fopen(specfile, 'rt');
if (specFID == -1)
    error(sprintf('Could not open spec file %s.', specfile));
end


line = fgetl(specFID);
dimensions = sscanf(line, '%i %i %i');
fclose(specFID);


