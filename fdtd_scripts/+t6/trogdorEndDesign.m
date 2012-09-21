function designObject = trogdorEndDesign(varargin)

simObject = t6.simulation().instance();

t6.TrogdorSimulation.clear();

designObject = t6.adjoint.DesignObject(simObject);
