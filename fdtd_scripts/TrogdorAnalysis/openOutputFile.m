% openOutputFile        Return FID and frame dimensions for FDTD output.
%   [fid, dimensions] = openOutputFile(fileprefix) reads frame dimensions
%   from the spec file fileprefix.txt and returns a handle to the binary
%   data from fileprefix.dat.
%
%   version 5.0
%   July 25, 2009

function [fid, dimensions] = openOutputFile(fileprefix, varargin)

dataFormat = 'ieee-le';
if (nargin > 1)
    dataFormat = varargin{1};
end

dimensions = [-1, -1, -1];
fid = -1;

datFile = [fileprefix, '.dat'];
specFile = [fileprefix, '.txt'];


specFID = fopen(specFile, 'rt');
if (specFID == -1)
    error(sprintf('Could not open spec file %s.', specFile));
end

line = fgetl(specFID);

% Trogdor 5 support is scanty here.
% A future solution will provide a more thorough way to access all kinds
% of data; the spec file should eventually also report endianness
% and precision.
if strcmp(line, 'trogdor5output')
    foundDimensions = 0;
    foundFields = 0;
    while ischar(line)
        line = fgetl(specFID);
        
        if ~ischar(line)
            break
        end
        
        tryDatafile = sscanf(line, 'datafile %s');
        if ~strcmp(tryDatafile, '')
            datFile = tryDatafile;
        end
        
        tryField = sscanf(line, 'field %s');
        if ~strcmp(tryField, '')
            foundFields = foundFields+1;
        end
        
        tryRegion = sscanf(line, ...
            'region [(%i, %i, %i), (%i, %i, %i)] stride (%i, %i, %i)');
        if length(tryRegion) == 9
            region = tryRegion;
            foundDimensions = foundDimensions+1;
        end
    end
    
    if foundDimensions == 0
        error(sprintf('Could not find dimensions of data file %s',...
            datFile));
    end
    
    if foundDimensions > 1
        error(sprintf('Output with more than one region is unsupported.'));
    end
    
    if length(region) ~= 9
        error(sprintf('Could not read region from %s', specFile));
    else
        dimensions = floor(region(4:6)./ region(7:9)) -...
            floor(region(1:3) ./ region(7:9)) + 1;
        dimensions(end+1) = foundFields;
    end
else
    dimensions = sscanf(line, '%i %i %i');
    datFile = [fileprefix, '.dat'];
end
fclose(specFID);

fid = fopen(datFile, 'r', 'ieee-le');
if (fid == -1)
    error(sprintf('Could not open data file %s.', datFile));
end

dimensions = transpose(dimensions);