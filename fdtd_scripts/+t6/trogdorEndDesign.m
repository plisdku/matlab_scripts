function designObject = trogdorEndDesign(varargin)

X.Tag = '';
X = parseargs(X, varargin{:});

simObject = t6.simulation();

designObject = t6.adjoint.DesignObject(simObject, X.Tag);

TROGDOR_SIMULATION = [];

%t6.TrogdorSimulation.clear();

