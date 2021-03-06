function xyzPos = positions(obj, varargin)
% positions  Get the positions of each sample in the output file
%
% xyzPos = file.positions returns a cell array of size [3, NumFields].
%
% Named parameters:
%
% Regions           'Separate' or 'Together'
% Region            index of region to get positions for
% Field             index of field to get positions for
% Size              3-element array [Nx Ny Nz], number of samples to take
% Bounds            Ask Paul
% InterpolateSpace  true or false
%

X.Regions = 'Separate'; % Separate or Together
X.Region = [];
X.Field = [];
X.Size = [];
X.Bounds = [];
X.InterpolateSpace = [];
X = parseargs(X, varargin{:});

if isempty(X.InterpolateSpace)
    X.InterpolateSpace = obj.hasBounds();
end

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
    xyzPos(:) = naturalSamplingPositions(obj, X.Bounds, X.Region, X.Size);
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
