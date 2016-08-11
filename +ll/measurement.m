function measurement(varargin)

X.Bounds = [];
X.F = [];
X.Jx = '0.0';
X.Jy = '0.0';
X.Jz = '0.0';
X.g = []; % can be a string or a callback
X.Filename = '';

X = parseargs(X, varargin{:});

extents = X.Bounds(4:6) - X.Bounds(1:3);

measStruct = struct('bounds', X.Bounds, ...
    'dimensions', nnz(extents), ...
    'F', X.F, 'Jx', X.Jx, 'Jy', X.Jy, 'Jz', X.Jz, 'g', X.g, ...
    'filename', X.Filename);

global LL_MODEL;
LL_MODEL.measurements{numel(LL_MODEL.measurements)+1} = measStruct;

