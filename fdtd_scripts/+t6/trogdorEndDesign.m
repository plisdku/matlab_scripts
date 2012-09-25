function designObject = trogdorEndDesign(varargin)

simObject = t6.simulation();

designObject = t6.adjoint.DesignObject(simObject);

TROGDOR_SIMULATION = [];

%t6.TrogdorSimulation.clear();

