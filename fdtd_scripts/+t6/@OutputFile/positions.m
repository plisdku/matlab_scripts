function xyzPos = positions(obj, varargin)

X.Regions = 'Separate'; % Separate or Together
X.Region = [];
X.Field = [];
X.Size = [];
X.Bounds = [];
X.InterpolateSpace = obj.hasBounds();
X = parseargs(X, varargin{:});

validateArguments(obj, X);
xyzPos = cell(3,1);

if isempty(X.Region), X.Region = 1; end
if isempty(X.Field), X.Field = 1; end

if strcmpi(X.Regions, 'Together') || ~X.InterpolateSpace
    offset = obj.Fields{X.Field}.Offset;
    yees = obj.yeeCells('Region', X.Region, 'Regions', X.Regions);
    
    for xyz = 1:3
        xyzPos{xyz} = obj.Dxyz(xyz)*(yees{xyz} + offset(xyz)) + obj.Origin(xyz);
    end
else
    xyzPos = naturalSamplingPositions(obj, X.Bounds, X.Region, X.Size);
end



function validateArguments(obj, X)

if ~all(islogical(X.InterpolateSpace))
    error('InterpolateSpace must be logical true or false');
end

if ~strcmpi(X.Regions, 'Separate') && ~strcmpi(X.Regions, 'Together')
    error('Regions must be Separate or Together');
end

if strcmpi(X.Regions, 'Separate')
    if obj.numRegions() > 1 && isempty(X.Region)
        error('More than one region is present in this file.  Please provide a Region argument.');
    end
    
    if X.InterpolateSpace && ~obj.hasBounds() && isempty(X.Bounds)
        error('To interpolate these fields please specify Bounds.');
    end
    
elseif strcmpi(X.Regions, 'Together')
    if ~isempty(X.Bounds)
        error('Bounds argument cannot be used with Regions = Together');
    end
    
    if ~isempty(X.Size)
        error('Size argument cannot be used with Regions = Together');
    end
end

if obj.numFields() > 1 && isempty(X.Field) && ~X.InterpolateSpace && ...
        strcmpi(X.Regions, 'Separate')
    error('More than one field  is present in this file.  Please provide a Field argument.');
end


function pos = naturalSamplingPositions(obj, bounds, region, numSamples)

pos = cell(1,3);

if isempty(bounds)
    bounds = obj.Regions.Bounds(region,:);
end

if isempty(numSamples)
    numSamples = (bounds(4:6) - bounds(1:3)) ./ obj.Dxyz;
end

numSamples = round(numSamples);
numSamples(numSamples < 1) = 1;

for xyz = 1:3
    if numSamples(xyz) > 1
        pos{xyz} = linspace(bounds(xyz), bounds(xyz+3), numSamples(xyz));
    else
        pos{xyz} = bounds(xyz);
    end
end    





%{
xx = cell(obj.numRegions(), length(obj.Fields));
yy = xx;
zz = xx;

[ii jj kk] = obj.yeeCells();

if ~iscell(ii)
    ii = {ii};
    jj = {jj};
    kk = {kk};
end

for rr = 1:obj.numRegions()
for ff = 1:length(obj.Fields)
    offset = obj.Fields{ff}.Offset;
    xx{rr,ff} = obj.Dxyz(1)*(ii{rr} + offset(1)) + obj.Origin(1);
    yy{rr,ff} = obj.Dxyz(2)*(jj{rr} + offset(2)) + obj.Origin(2);
    zz{rr,ff} = obj.Dxyz(3)*(kk{rr} + offset(3)) + obj.Origin(3);
end
end

if length(xx) == 1
    xx = xx{1};
    yy = yy{1};
    zz = zz{1};
end
%}