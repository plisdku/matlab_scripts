function outStruct = solveTE(varargin)
% propagation along z
% parallel wavevector along y

import tmm.*

X.boundaryZ = [];
X.epsr = [];
X.mur = [];
X.omega = [];
X.ky = [];
X.sourceExLeft = 0;
X.sourceExRight = 0;
X.Jx = []; % Source currents at each boundary position
X.My = []; % Source currents at each boundary position
X.Ex = []; % Output positions for Ex etc.
X.Exy = []; % dEx/dy positions
X.Exz = []; % dEx/dz positions
X.Hy = [];
X.Hyy = [];
X.Hyz = [];
X.Hz = [];
X.Hzy = [];
X.Hzz = [];
X.forceBoundModes = false;
X = parseargs(X, varargin{:});

if isempty(X.Jx)
    X.Jx = zeros(1, numel(X.boundaryZ));
end

if isempty(X.My)
    X.My = zeros(1, numel(X.boundaryZ));
end

outStruct = solveTE_fast(X.boundaryZ, X.epsr, X.mur, X.omega, X.ky, ...
    X.sourceExLeft, X.sourceExRight, ...
    X.Jx, X.My, ...
    X.Ex, X.Exy, X.Exz, ...
    X.Hy, X.Hyy, X.Hyz, ...
    X.Hz, X.Hzy, X.Hzz, ...
    X.forceBoundModes);

