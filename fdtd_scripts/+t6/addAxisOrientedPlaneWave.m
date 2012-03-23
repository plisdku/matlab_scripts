function addAxisOrientedPlaneWave(varargin)
%addAxisOrientedPlaneWave Source fields into the grid with a plane wave
%   addAxisOrientedPlaneWave('YeeCells', [0 0 0 100 100 100], ...
%       'Field', 'ex', 'Direction', [0 0 1], 'TimeData', sin(1:numT)) will add
%       an ex-polarized plane wave with time dependence sin(timestep) into the
%       grid with a total-field scattered-field boundary.
%
%   Usage: addAxisOrientedPlaneWave(named parameters)
%
%   Named parameters:
%       YeeCells    The cells making up the total-field region of the grid;
%                   [x0 y0 z0 x1 y1 z1] specifies all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.
%                   (YeeCells or Bounds required)
%       Bounds      The real extent of the total-field region of the grid;
%                   [x0 y0 z0 x1 y1 z1] will be rounded outwards to the next
%                   largest full Yee rectangle [m0 n0 p0 m1 n1 p1], suitably for
%                   the grid resolution.
%                   (YeeCells or Bounds required)
%       Duration    The range of timesteps on which to source fields; [t0 t1]
%                   will source on timesteps t such that t0 <= t <= t1.  Using
%                   multiple rows specifies multiple ranges of timesteps.
%                   (default: all timesteps)
%       FieldFunction   A function of time, e.g. @(t) sin(omega*t).  If more
%                   than one field component is specified, then this argument
%                   must be a cell array with one function per field component,
%                   e.g. {exFunc eyFunc}.
%                   (FieldFunction or TimeData required)
%       TimeData    An array of size [nFields nTimesteps].  If the Duration
%                   is specified as [0 10] then TimeData needs 11 columns, one
%                   for each sourced timestep.  (required)
%       OmitSide    An axis-aligned unit vector or cell array of axis-aligned
%                   unit vectors.  The unit vector [1 0 0] specifies the +X
%                   side of the TFSF boundary and instructs Trogdor to not
%                   add or subtract electromagnetic fields on that side.
%                   (default: empty cell array)
%
%   Example:
%
%   addAxisOrientedPlaneWave('YeeCells', [0 0 0 100 100 100], 'TimeData', ...
%       cos(1:numT), 'Field', 'ex', 'Direction', [0 0 1], ...
%       'OmitSide', {[0 0 1], [0 -1 0]});
%
%   Source data into the box [0 0 0 100 100 100] but omit the TFSF correction on
%   the +Z and -Y faces.

import t6.*

grid = t6.TrogdorSimulation.instance().currentGrid();

X.Field = '';
X.YeeCells = [];
X.Bounds = [];
X.Duration = [];
X.Direction = [];
X.FieldFunction = [];
X.TimeData = [];
X.OmitSide = [];
X = parseargs(X, varargin{:});

t6.validateYeeCellsAndBounds(X);

% Validate fields; should be a single string with some tokens in it
fieldTokens = tokenizeFields(X.Field, 'e h');

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
if ~isempty(X.Bounds)
    X.YeeCells = boundsToYee(X.Bounds, fieldTokens);
end

if length(X.Duration) == 0
    X.Duration = [0, t6.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

numT = sum(X.Duration(:,2) - X.Duration(:,1) + 1);

%  Validate time data.
if length(X.TimeData) ~= 0
    if length(X.TimeData) ~= numT
        error('TimeData must have the same length as the Duration');
    elseif size(X.TimeData, 1) ~= length(fieldTokens) && ...
        size(X.TimeData, 1) ~= 1
        error('TimeData must have size [nFields, timesteps] or [timsteps].');
    end
end

% Validate direction
if length(X.Direction) ~= 3
    error('Direction must be an axis-oriented unit vector, e.g. [1 0 0].');
end
if norm(X.Direction) ~= 1
    error('Direction must be an axis-oriented unit vector.');
end
if abs(X.Direction(1)) ~= 1 && abs(X.Direction(2)) ~= 1 && ...
    abs(X.Direction(3)) ~= 1
    error('Direction must be an axis-oriented unit vector.');
end

obj = struct;
obj.type = 'TFSFSource';
obj.field = fieldTokens;
obj.yeeCells = X.YeeCells;
obj.duration = X.Duration;
obj.timeData = X.TimeData;
obj.fieldFunction = X.FieldFunction;
obj.direction = X.Direction;

% Validate omitted sides, kind of
if isempty(X.OmitSide)
    obj.omitSides = {};
elseif ~iscell(X.OmitSide)
    obj.omitSides = {X.OmitSide};
else
    obj.omitSides = X.OmitSide;
end


grid.TFSFSources = {grid.TFSFSources{:}, obj};
