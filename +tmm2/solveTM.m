function outStruct = solveTM(varargin)
% propagation along z
% parallel wavevector along y

import tmm.*

X.boundaryZ = [];
X.epsr = [];
X.mur = [];
X.omega = [];
X.ky = [];
X.sourceHxLeft = 0;
X.sourceHxRight = 0;
X.Mx = []; % Source currents at each boundary position
X.Jy = []; % Source currents at each boundary position
X.Hx = []; % Output positions for Hx etc.
X.Hxy = []; % dHx/dy positions
X.Hxz = []; % dHx/dz positions
X.Ey = [];
X.Eyy = [];
X.Eyz = [];
X.Ez = [];
X.Ezy = [];
X.Ezz = [];
X.forceBoundModes = false;
X = parseargs(X, varargin{:});

if isempty(X.Mx)
    X.Mx = zeros(1, numel(X.boundaryZ));
end

if isempty(X.Jy)
    X.Jy = zeros(1, numel(X.boundaryZ));
end

outStruct = solveTM_fast(X.boundaryZ, X.epsr, X.mur, X.omega, X.ky, ...
    X.sourceHxLeft, X.sourceHxRight, ...
    X.Mx, X.Jy, ...
    X.Hx, X.Hxy, X.Hxz, ...
    X.Ey, X.Eyy, X.Eyz, ...
    X.Ez, X.Ezy, X.Ezz, ...
    X.forceBoundModes);

