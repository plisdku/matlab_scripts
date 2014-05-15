function designObject = trogdorEndDesign(varargin)
% designObject = trogdorEndDesign(varargin)
%
% Tag     simulation tag for output and sim directories

X.Tag = '';
X = parseargs(X, varargin{:});

simObject = t6.simulation();

designObject = t6.adjoint.DesignObject(simObject, X.Tag);

TROGDOR_SIMULATION = [];

%t6.TrogdorSimulation.clear();

