function data = readFrames(obj, varargin)
% outputFile.readFrames()
% outputFile.readFrames(numFrames)

X.NumFrames = 1;
X.Regions = 'Separate';
X = parseargs(X, varargin{:});

if obj.FileHandle == -1
    error('No data file is open.  Try open()?');
end

if nargin > 1
    numFrames = X.NumFrames;
else
    numFrames = 1;
end

if strcmp(X.Regions, 'Separate')
    data = readFrames_SeparateRegions(obj, numFrames);
else
    data = readFrames_RegionsTogether(obj, numFrames);
end

