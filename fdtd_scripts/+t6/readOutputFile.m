function [data, positions] = readOutputFile(fileName, varargin)
%readOutputFile Grab all the data from the named output file
%   alldat = readOutputFile(filename) is a shortcut for using OutputFile to
%   read the file.
%
% Regions
% Positions
% Size
% Times
% InterpolateSpace

X.Regions = 'Separate';
X.Positions = [];
X.Size = [];
X.Times = [];
X.InterpolateSpace = [];
X = parseargs(X, varargin{:});

file = t6.OutputFile(fileName);
%try
file.open();

if nargout > 1
    [data, positions] = file.readFrames(...
        'Regions', X.Regions, ...
        'Positions', X.Positions, ...
        'Times', X.Times, ...
        'Size', X.Size, ...
        'InterpolateSpace', X.InterpolateSpace);
else
    data = file.readFrames(...
        'Regions', X.Regions, ...
        'Positions', X.Positions, ...
        'Times', X.Times, ...
        'Size', X.Size, ...
        'InterpolateSpace', X.InterpolateSpace);
end
    
%end
file.close();



