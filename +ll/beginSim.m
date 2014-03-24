function beginSim(varargin)

X.Bounds = [];
X.Frequency = [];
X.PML = [];
X.HMin = '';
X.HMax = '';

X = parseargs(X , varargin{:});

global LL_MODEL;

LL_MODEL = struct('bounds', X.Bounds, ...
    'frequency', X.Frequency, ...
    'PML', X.PML, ...
    'PMLBounds', X.Bounds - [1 1 1 -1 -1 -1].*X.PML, ...
    'PMLThickness', max(abs(X.PML)), ...
    'sources', {{}}, ...
    'outputs', {{}}, ...
    'measurements', {{}}, ...
    'materials', {{}}, ...
    'incidentField', struct('Ex', '', 'Ey', '', 'Ez', ''), ...
    'meshes', {{}}, ...
    'hmin', X.HMin, ...
    'hmax', X.HMax);
