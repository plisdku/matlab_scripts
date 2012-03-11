function [data, positions] = readFrames(obj, varargin)
% outputFile.readFrames()
% outputFile.readFrames(numFrames)
% outputFile.readFrames(numFrames, 'Regions', 'Together')
% outputFile.readFrames(numFrames, 'Regions', 'Separate')

X.NumFrames = obj.numFramesAvailable();
X.Regions = 'Separate';
X.Size = [];
X.Times = [];
X.Positions = [];
X.InterpolateSpace = [];

X = parseargs(X, varargin{:});
if isempty(X.InterpolateSpace)
    X.InterpolateSpace = obj.hasBounds();
end

if obj.FileHandle == -1
    error('No data file is open.  Try open()?');
end

if ~isempty(X.Positions)
    if ~iscell(X.Positions)
        X.Positions = {X.Positions(:,1), X.Positions(:,2), X.Positions(:,3)};
    end
end

% Figure out whether to interpolate in space and/or time.
% By default, as much interpolation as possible will be performed.
% The user can turn it off of course.

if X.InterpolateSpace
    if isempty(X.Positions)
        % fill in some natural sampling points.
        X.Positions = myNaturalSamplingPositions(obj, X.Size);
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

if ~isempty(X.Positions)
    positions = X.Positions;
elseif nargout > 1
    if obj.numFields() == 1
        positions = obj.positions();
    else
        warning('There are multiple fields in this file.  Not computing positions.');
    end
end




function pos = myNaturalSamplingPositions(obj, numSamplesInRegions)

pos = cell(obj.numRegions, 3);

for rr = 1:obj.numRegions()
    
    if ~isempty(numSamplesInRegions)
        numSamples = numSamplesInRegions(rr,:);
    else
        numSamples = [];
    end
    
    positionsInRegion = naturalSamplingPositions(obj, obj.Regions.Bounds(rr,:), rr, ...
        numSamples);
    
    for xyz = 1:3
        pos{rr,xyz} = positionsInRegion{xyz};
    end
end


%{
pos = cell(obj.numRegions, 3);
for rr = 1:obj.numRegions()
    if isempty(numSamplesInRegions)
        numSamples = (obj.Regions.Bounds(rr,4:6) - obj.Regions.Bounds(rr,1:3)) ...
            ./ obj.Dxyz;
    else
        numSamples = numSamplesInRegions(rr,:);
    end
    numSamples = ceil(numSamples);
    numSamples(numSamples < 1) = 1;

    for xyz = 1:3
        if numSamples(xyz) > 1
            pos{rr,xyz} = linspace(obj.Regions.Bounds(rr,xyz), ...
                obj.Regions.Bounds(rr,xyz+3), numSamples(xyz));
        else
            pos{rr,xyz} = obj.Regions.Bounds(rr,xyz);
        end
    end    
end
%}


