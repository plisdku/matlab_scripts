function outStruct = solveTE(varargin)
% propagation along z
% parallel wavevector along y

import tmm2.*

X.boundaryX = [];
X.epsr = [];
X.mur = [];
X.omega = [];
X.ky = [];
X.sourceEzLeft = 0;
X.sourceEzRight = 0;
X.Jz = []; % Source currents at each boundary position
X.Mx = []; % Source currents at each boundary position
X.Ez = []; % Output positions for Ex etc.
X.Ezx = []; % dEx/dy positions
X.Ezy = []; % dEx/dz positions
X.Hx = [];
X.Hxx = [];
X.Hxy = [];
X.Hy = [];
X.Hyx = [];
X.Hyy = [];
X.forceBoundModes = false;
X = parseargs(X, varargin{:});

if isempty(X.Jz)
    X.Jz = zeros(1, numel(X.boundaryX));
end

if isempty(X.Mx)
    X.Mx = zeros(1, numel(X.boundaryX));
end

outStruct = solveTE_fast(X.boundaryX, X.epsr, X.mur, X.omega, X.ky, ...
    X.sourceEzLeft, X.sourceEzRight, ...
    X.Jz, X.Mx, ...
    X.Ez, X.Ezx, X.Ezy, ...
    X.Hx, X.Hxx, X.Hxy, ...
    X.Hy, X.Hyx, X.Hyy, ...
    X.forceBoundModes);

