function measurement(varargin)

X.Bounds = [];
X.F = [];
X.Jx = '0.0';
X.Jy = '0.0';
X.Jz = '0.0';
X.g = '0.0';

X = parseargs(X, varargin{:});

measStruct = struct('bounds', X.Bounds, ...
    'F', X.F, 'Jx', X.Jx, 'Jy', X.Jy, 'Jz', X.Jz, 'g', X.g);

global LL_MODEL;
LL_MODEL.measurements{numel(LL_MODEL.measurements)+1} = measStruct;

