function incidentField(varargin)

X.Ex = '';
X.Ey = '';
X.Ez = '';

X = parseargs(X, varargin{:});

incidentFieldStruct = struct('Ex', X.Ex, 'Ey', X.Ey, 'Ez', X.Ez);

global LL_MODEL;
LL_MODEL.incidentField = incidentFieldStruct;