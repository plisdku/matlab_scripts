function addHardSource(varargin)
%addHardSource Add a hard electromagnetic field source to the current grid.
%   addHardSource('Field', 'ez', 'YeeCells', [50 50 50 50 50 50], 'TimeData', ...
%       sin(1:numT)) will set the ez field at (50,50,50) to sin(timestep) on
%       every timestep.
%
%   Usage: addHardSource(named parameters)
%
%   Named parameters:
%       Field       The electromagnetic fields to source.  Valid fields:
%                   'ex', 'ey', 'ez'    electric fields
%                   'hx', 'hy', 'hz'    magnetic fields
%                   Any combination of fields may be used, but they should be
%                   specified in order: ex before ey, electric before magnetic.
%                   (required)
%       YeeCells    The region of the grid in which to add electromagnetic
%                   field; [x0 y0 z0 x1 y1 z1] will source all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.  Multiple
%                   rows may be used to source in multiple regions.
%                   (required)
%       Duration    The range of timesteps on which to source fields; [t0 t1]
%                   will source on timesteps t such that t0 <= t <= t1.  Using
%                   multiple rows specifies multiple ranges of timesteps.
%                   (default: [0, numT - 1] for a simulation with numT timesteps)
%       TimeData    Fields on each timestep, an array of size [nFields nTimesteps].
%                   If the Duration is specified as [0 10] then TimeData needs 11
%                   columns, one for each sourced timestep.
%                   Units: SI (e in V/m, h in A/m).
%                   (TimeData or SpaceTimeData required)
%       MaskData    A per-field, per-cell prefactor to multiply the source by.
%                   If YeeCells specifies a rectangle of size [nx ny nz] then
%                   MaskData must have size [nx ny nz nFields].  If YeeCells
%                   has multiple rows, then MaskData must be a cell array 
%                   where each cell contains the mask for one row in YeeCells.
%                   (default: no mask)
%       SpaceTimeData   Specifies the field for each cell at each timestep.  If
%                   YeeCells specifies a rectangle of size [nx ny nz] then
%                   SpaceTimeData must have size [nx ny nz nFields nTimesteps].
%                   If YeeCells has multiple rows, then SpaceTimeData must be
%                   a cell array with one entry per row of YeeCells.  If the 
%                   Duration is specified as [0 10] then SpaceTimeData must
%                   provide data for 11 timesteps.  See TimeData for units.
%                   (TimeData or SpaceTimeData required.)
%
%   Example:
%
%   addHardSource('Field', 'ex hx', 'YeeCells', [50 50 50 50 50 50; ...
%       60 60 60 60 60 60], 'Duration', [0 10; 100 110], 'TimeData', ...
%       [sin(0:10), cos(100:110); exp(-(0:10)), exp(-(100:110))])
%
%   Specify soft source in cells (50, 50, 50) and (60, 60, 60), nonzero for
%   timesteps 0:10 and 100:110, where ex varies sinusoidally and hx varies
%   as a decaying exponential.
grid = t5.TrogdorSimulation.instance().currentGrid();

X.Field = '';
X.YeeCells = [];
X.Duration = [];
X.TimeData = [];
X.MaskData = [];
X.SpaceTimeData = [];
X = parseargs(X, varargin{:});

t5.validateSourceDataParameters(X); % will call error() for problems

% Now I need to validate the Field: it must be ex, ey, ez, hx, hy, hz stuff.
fieldTokens = {};
remainder = X.Field;
while ~strcmp(remainder, '')
    [token, remainder] = strtok(remainder);
    if ~strcmp(token, '')
        fieldTokens = {fieldTokens{:}, token};
        
        if length(token) ~= 2
            error('Bad field %s', token);
        elseif (token(1) ~= 'e' && token(1) ~= 'h') || ...
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

% Validate mask
if length(X.MaskData) ~= 0
    if ~iscell(X.MaskData)
        X.MaskData = {X.MaskData};
    end
    if length(X.MaskData) ~= size(X.YeeCells, 1)
        error('Please provide each source region mask in its own cell.');
    end
    
    for ll = 1:length(X.MaskData)
        yeeSize = X.YeeCells(ll,4:6) - X.YeeCells(ll,1:3) + 1;
        if size(X.MaskData{ll}) ~= [yeeSize, length(fieldTokens)]
            error('Mask #%i needs to have size [yeeCells, numFields]', ll);
        end
    end
end

% Validate space-time data
if length(X.SpaceTimeData) ~= 0
    if ~iscell(X.SpaceTimeData)
        X.SpaceTimeData = {X.SpaceTimeData};
    end
    if length(X.SpaceTimeData) ~= size(X.YeeCells, 1)
        error('Please provide each source region spacetime data in its own cell.');
    end
    
    for ll = 1:length(X.SpaceTimeData)
        yeeSize = X.YeeCells(ll,4:6) - X.YeeCells(ll,1:3) + 1;
        if size(X.SpaceTimeData{ll}) ~= [yeeSize, length(fieldTokens), numT]
            error('Space time data cell #%i is the wrong size.', ll);
        end
    end
end


obj = struct;
obj.type = 'HardSource';
obj.field = fieldTokens;
obj.yeeCells = X.YeeCells;
obj.duration = X.Duration; % validated
obj.timeData = X.TimeData; % validated
obj.maskData = X.MaskData; % validated
obj.spaceTimeData = X.SpaceTimeData;

grid.HardSources = {grid.HardSources{:}, obj};
