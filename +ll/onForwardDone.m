function onForwardDone(varargin)

X.Callback = [];

X = parseargs(X, varargin{:});

if ~isa(X.Callback, 'function_handle')
    error('Callback must be a function handle');
end

global LL_MODEL;

if ~isequal(LL_MODEL.forwardCallback, [])
    warning('onForwardDone callback is being overwritten');
end

LL_MODEL.forwardCallback = X.Callback;


