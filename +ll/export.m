function export(varargin)

X.Mode = {'Forward', 'Adjoint'};
X.File = '';
X.Expression = '';
X.Points = [];

X = parseargs(X, varargin{:});

global LL_MODEL;

exportStruct = struct('mode', X.Mode, 'file', X.File, ...
    'expression', X.Expression, 'points', X.Points);

LL_MODEL.exports{numel(LL_MODEL.exports)+1} = exportStruct;

