function setPML(varargin)
%setPML Adjust the depth and variable parameters of the CFS-RIPML for a grid or
%   for a particular material
%
%   setPML('Depth', [10 10 0 10 10 0]) will add 10 cells of PML to the current
%       grid on the -X, -Y, +X and +Y sides
%   setPML('Sigma', '(d^4)*0.8*4/(((mu0/eps0)^0.5)*dx') will cause the PML to
%       use a fourth-power scaling for the conductivity, where the variable 'd'
%       varies from 0 at the inside edge of the PML to 1 at the deepest point,
%       and 'dx' takes on the value of dx for +X and -X PML, dy for +Y and -Y,
%       and dz for +Z and -Z.
%   setPML('Kappa', '5') will cause the PML to use a fixed real stretch factor
%       of 5.
%   setPML('Alpha', '0') will turn off the CFS-RIPML's complex frequency shift.
%   setPML('Material', 'Gold', 'Alpha', '0') will set alpha to zero for gold
%       PML on all grids.
%
%   Named parameters may be combined, so 'Sigma', 'Alpha', 'Kappa' and 'Depth'
%   may be set all at once; however, 'Depth' and 'Material' cannot be set at
%   once as all materials have the same PML thickness.

sim = t5.TrogdorSimulation.instance();
grid = t5.TrogdorSimulation.instance().currentGrid();

X.Depth = [];
X.Material = [];
X.Kappa = '';
X.Alpha = '';
X.Sigma = '';
X = parseargs(X, varargin{:});

if length(X.Depth) ~= 0
    if length(X.Depth) ~= 6
        error('PML depth must be a length-six array of integers.');
    end
    
    if length(X.Material) ~= 0
        error('Cannot set PML depth per material.');
    end
    grid.PML = X.Depth;
end

if any(X.Depth < 0)
    error('PML depths must all be nonnegative.');
end

if length(X.Material) ~= 0
    
    index = t5.indexOf(X.Material, sim.Materials);
    if index == -1
        error('Material not found.');
    end
    
    if length(X.Kappa) ~= 0
        sim.Materials{index}.PMLParams.kappa = X.Kappa;
    end
    
    if length(X.Sigma) ~= 0
        sim.Materials{index}.PMLParams.sigma = X.Sigma;
    end
    
    if length(X.Alpha) ~= 0
        sim.Materials{index}.PMLParams.alpha = X.Alpha;
    end
else
    if length(X.Kappa) ~= 0
        grid.PMLParams.kappa = X.Kappa;
    end
    
    if length(X.Sigma) ~= 0
        grid.PMLParams.sigma = X.Sigma;
    end
    
    if length(X.Alpha) ~= 0
        grid.PMLParams.alpha = X.Alpha;
    end
end
