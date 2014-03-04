function output(varargin)

X.Bounds = [];

X = parseargs(X, varargin{:});

outputStruct = struct('bounds', X.Bounds);

global LL_MODEL;
LL_MODEL.outputs{numel(LL_MODEL.outputs)+1} = outputStruct;