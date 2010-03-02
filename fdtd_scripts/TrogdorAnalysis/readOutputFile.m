% readOutputFile        Grab all the data from the named output file
%   alldat = readOutputFile(filename) is a shortcut for using
%   openOutputFile, getOutputFrames and fclose.
%
%   See also: openOutputFile, getOutputFrames, outputFFT, outputHarmonic
%
%   version 4.5
%   July 29, 2008

function alldat = readOutputFile(filename)

[fid, dim] = openOutputFile(filename);
alldat = getOutputFrames(fid, dim);
fclose(fid);
