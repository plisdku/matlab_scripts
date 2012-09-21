function addCurrentSource(varargin)
%addCurrentSource Add a source of electric and/or magnetic current density
% to the grid.
%
% addCurrentSource('Field', 'my', 'Bounds', [0 0 0 0 100 0], ...
%   'FieldFunction', @(x,y,z,t) sin(2*pi*y/200)*sin(t/10))
% Drive the grid with a sinusoidal magnetic current density along a line
%   from (0,0,0) to (0,100,0).
%
%   Usage: addCurrentSource(named parameters)
%
%   Named parameters:
%       Field       The electromagnetic fields to source.  Valid fields:
%                   'mx', 'my', 'mz'    magnetic current densities
%                   'jx', 'jy', 'jz'    electric current densities
%                   Any combination of fields may be used, but they should be
%                   specified in order: mx before my, electric before
%                   magnetic.
%                   (required)
%       YeeCells    The region of the grid in which to add electromagnetic
%                   current; [x0 y0 z0 x1 y1 z1] will source all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.  Multiple
%                   rows may be used to source in multiple regions.
%                   (YeeCells or Bounds required)
%       Bounds      The region of the simulation space in which to add
%                   electromagnetic current; [x0 y0 z0 x1 y1 z1] in real units
%                   will be used to choose the YeeCells [m0 n0 p0 m1 n1 p1] in
%                   which to source fields, suitably for the grid resolution
%                   (YeeCells or Bounds required)
%       Duration    The range of timesteps on which to source currents; [t0 t1]
%                   will source on timesteps t such that t0 <= t <= t1.  Using
%                   multiple rows specifies multiple ranges of timesteps.
%                   (default: all timesteps)
%       FieldFunction   A function of space and time, e.g.
%                   @(x, y, z, t) sin(x).*cos(y).*tan(z).*exp(-t).
%                   If more than one current component is specified, then this
%                   argument must be a cell array with one function per field
%                   component, e.g. {jxFunc, jyFunc, jzFunc}.
%                   (FieldFunction, TimeData, or SpaceTimeData required)
%       FieldFunctor    A function F(x,y,z) that returns a function f(t)
%                   providing the field on each timestep.
%       TimeData    An array of size [nFields nTimesteps].  If the Duration
%                   is specified as [0 10] then TimeData needs 11 columns, one
%                   for each sourced timestep.
%                   (FieldFunction, TimeData or SpaceTimeFile required)
%       SpaceTimeData An array of size [nFields nTimesteps*nCells].
%       SpaceTimeFile   Name of a file that will contain the current for each
%                   cell at each timestep.  When using SpaceTimeFile, Trogdor
%                   will run once to generate a data request file named
%                   [SpaceTimeFile, '.m'] containing an ordered list of points
%                   and a set of timesteps at which to provide electromagnetic
%                   current values.  The space-time data file must be provided
%                   according to this data request before the simulation may
%                   run.  (SpaceTimeFile or TimeData required.)
%
%   Example:
%
%   addCurrentSource('Field', 'jz', 'Bounds', [50 50 50 50 50 50], ...
%       'FieldFunction', @(x,y,z,t) sin(t/20));
%
%   Specify a point z-oriented dipole of amplitude 1 at (50,50,50).

import t6.*

sim = simulation();

X.Field = '';
X.YeeCells = [];
X.Bounds = [];
X.Duration = [];
X.FieldFunction = [];
X.FieldFunctor = [];
X.TimeData = [];
X.SpaceTimeData = [];
X.Mode = '';
%X.SpaceTimeFile = [];
X = parseargs(X, varargin{:});

validateDataRequestParameters(X);
validateYeeCellsAndBounds(X);

% Validate fields; should be a single string with some tokens in it
fieldTokens = tokenizeFields(X.Field, 'j m je mh');

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
if ~isempty(X.Bounds)
    X.YeeCells = sim.boundsToYee(X.Bounds, fieldTokens);
end

% Validate duration
if isempty(X.Duration)
    X.Duration = [0, sim.NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

% Evaluate the field function at all the right places and times
if ~isempty(X.FieldFunction)
    if ~iscell(X.FieldFunction)
        X.FieldFunction = {X.FieldFunction};
    end
end

if ~isempty(X.FieldFunctor)
    if ~iscell(X.FieldFunctor)
        X.FieldFunctor = {X.FieldFunctor};
    end
end

if ~isempty(X.Mode)
    if ~strcmpi(X.Mode, 'forward') && ~strcmpi(X.Mode, 'adjoint')
        error('Mode must be forward or adjoint if given');
    end
end


obj = struct;
obj.type = 'CurrentSource';
obj.field = fieldTokens;
obj.yeeCells = X.YeeCells;
obj.bounds = X.Bounds;
obj.duration = X.Duration;
obj.fieldFunction = X.FieldFunction;
obj.fieldFunctor = X.FieldFunctor;
obj.timeData = X.TimeData;
obj.spaceTimeData = X.SpaceTimeData;
obj.mode = X.Mode;

%obj.spaceTimeFile = X.SpaceTimeFile;

sim.CurrentGrid.CurrentSources{end+1} = obj;


