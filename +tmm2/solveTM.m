function outStruct = solveTM(varargin)
% propagation along x
% parallel wavevector along y

import tmm2.*

X.boundaryX = [];
X.epsr = [];
X.mur = [];
X.omega = [];
X.ky = [];
X.sourceHzLeft = 0;
X.sourceHzRight = 0;
X.Mz = []; % Source currents at each boundary position
X.Jy = []; % Source currents at each boundary position
X.Hz = []; % Output positions for Hx etc.
X.Hzx = []; % dHx/dy positions
X.Hzy = []; % dHx/dz positions
X.Ex = [];
X.Exx = [];
X.Exy = [];
X.Ey = [];
X.Eyx = [];
X.Eyy = [];
X.Bz = [];
X.Dx = [];
X.Dy = [];
X.forceBoundModes = false;
X = parseargs(X, varargin{:});

if isempty(X.Mz)
    X.Mz = zeros(1, numel(X.boundaryX));
end

if isempty(X.Jy)
    X.Jy = zeros(1, numel(X.boundaryX));
end

outStruct = solveTM_fast(X.boundaryX, X.epsr, X.mur, X.omega, X.ky, ...
    X.sourceHzLeft, X.sourceHzRight, ...
    X.Mz, X.Jy, ...
    X.Hz, X.Hzx, X.Hzy, ...
    X.Ex, X.Exx, X.Exy, ...
    X.Ey, X.Eyx, X.Eyy, ...
    X.Bz, X.Dx, X.Dy, ...
    X.forceBoundModes);

