function addCurrentSource(varargin)
%addCurrentSource Add a source of electric and/or magnetic current density to the
%   grid.  addCurrentSource('Field', 'ky', 'YeeCells', [0 0 0 0 100 0], 'TimeData', ...
%       sin(1:numT)) will drive the grid with a sinusoidal magnetic current ky
%       along a line from (0,0,0) to (0,100,0) on every timestep.  The
%       units of current density are SI, i.e. J in amperes/m^2, K in V/m^2.
%
%   The electric and magnetic current densities are applied in the update equations as
%       dD/dt = curl(H) - J
%       dB/dt = -curl(E) - K
%
%   A single cell with a Jz source may be thought of as a dx*dy*dz rectangle of
%   constant current density.
%   
%   Usage: addCurrentSource(named parameters)
%
%   Named parameters:
%       Field       The electromagnetic fields to source.  Valid fields:
%                   'kx', 'ky', 'kz'    magnetic currents
%                   'jx', 'jy', 'jz'    electric currents
%                   Any combination of fields may be used, but they should be
%                   specified in order: kx before ky, electric before magnetic
%                   (required)
%       YeeCells    The region of the grid in which to add electromagnetic
%                   current; [x0 y0 z0 x1 y1 z1] will source all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.  Multiple
%                   rows may be used to source in multiple regions.
%                   (required)
%       Duration    The range of timesteps on which to source currents; [t0 t1]
%                   will source on timesteps t such that t0 <= t <= t1.  Using
%                   multiple rows specifies multiple ranges of timesteps.
%                   (default: [0, numT - 1] for a simulation with numT timesteps)
%       TimeData    Current density on each timestep, an array of size
%                   [nFields nTimesteps].  If the Duration is specified as [0 10]
%                   then TimeData needs 11 columns, one for each sourced timestep.
%                   Units: SI (j in A/m^2, k in V/m^2)
%                   (TimeData or SpaceTimeFile required)
%       MaskFile    Name of a file that will contain a per-field, per-cell
%                   prefactor to multiply the source by.  When using MaskFile,
%                   Trogdor will run once to generate a data request file named
%                   [MaskFile, '.m'] containing an ordered list of points at
%                   which to provide mask data.  The mask data file must be
%                   provided according to the data request before the simulation
%                   may run.  (default: no mask file)
%       SpaceTimeFile   Name of a file that will contain the current for each
%                   cell at each timestep.  When using SpaceTimeFile, Trogdor
%                   will run once to generate a data request file named
%                   [SpaceTimeFile, '.m'] containing an ordered list of points
%                   and a set of timesteps at which to provide electromagnetic
%                   current values.  The space-time data file must be provided
%                   according to this data request before the simulation may
%                   run.  See TimeData for units in the file.
%                   (SpaceTimeFile or TimeData required.)
%
%   Example:
%
%   addCurrentSource('Field', 'jz', 'YeeCells', [50 50 50 50 50 50; ...
%       60 60 60 60 60 60], 'Duration', [0 10; 100 110], 'SpaceTimeFile', 
%       'currentInput');
%
%   Specify a source of jz electric current in cells (50,50,50) and (60,60,60).
%   The first time Trogdor is run it will generate a data request file called
%   'currentInput.m' containing the ordered list of fields, cells and timesteps
%   at which data is required.  After this data is written in the correct order
%   to a binary file called 'currentInput', Trogdor can be run again to
%   completion.
grid = t5.TrogdorSimulation.instance().currentGrid();

X.Field = '';
X.YeeCells = [];
X.Duration = [];
X.TimeData = [];
X.MaskFile = [];
X.SpaceTimeFile = [];
X = parseargs(X, varargin{:});

t5.validateDataRequestParameters(X);

% Validate fields; should be a single string with some tokens in it
fieldTokens = {};
remainder = X.Field;
while ~strcmp(remainder, '')
    [token, remainder] = strtok(remainder);
    if ~strcmp(token, '')
        fieldTokens = {fieldTokens{:}, token};
        
        if length(token) ~= 2
            error('Bad field %s', token);
        elseif (token(1) ~= 'j' && token(1) ~= 'k') || ...
            (token(2) < 'x' || token(2) > 'z')
            token(1)
            token(2)
            error('Bad field %s', token);
        end
    end
end

% Validate duration
if length(X.Duration) == 0
    X.Duration = [0, t5.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

obj = struct;
obj.type = 'CurrentSource';
obj.field = fieldTokens;
obj.yeeCells = X.YeeCells;
obj.duration = X.Duration;
obj.timeData = X.TimeData;
obj.maskFile = X.MaskFile;
obj.spaceTimeFile = X.SpaceTimeFile;

grid.CurrentSources = {grid.CurrentSources{:}, obj};
