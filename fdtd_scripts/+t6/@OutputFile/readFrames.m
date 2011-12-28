function data = readFrames(obj, varargin)
% outputFile.readFrames()
% outputFile.readFrames(numFrames)
% outputFile.readFrames(numFrames, 'Regions', 'Together')
% outputFile.readFrames(numFrames, 'Regions', 'Separate')

X.NumFrames = obj.numFramesAvailable();
X.Regions = 'Separate';
X.Size = [];
X.Times = [];
X.Positions = [];
X.InterpolateSpace = obj.hasBounds();

X = parseargs(X, varargin{:});

if obj.FileHandle == -1
    error('No data file is open.  Try open()?');
end

% Figure out whether to interpolate in space and/or time.
% By default, as much interpolation as possible will be performed.
% The user can turn it off of course.

if X.InterpolateSpace
    if isempty(X.Positions)
        % fill in some natural sampling points
    end
elseif ~isempty(X.Positions)
    X.InterpolateSpace = true;
end

interpTime = ~isempty(X.Times);

if strcmp(X.Regions, 'Separate')
    if interpTime
        data = obj.readFrames_RegionsSeparate_InterpTime(X);
    else
        data = obj.readFrames_RegionsSeparate(X);
    end
else
    data = obj.readFrames_RegionsTogether(X.NumFrames);
end

