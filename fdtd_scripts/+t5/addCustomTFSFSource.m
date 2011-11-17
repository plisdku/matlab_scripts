function addCustomTFSFSource(varargin)
%addCustomTFSFSource Source fields into the grid with a data request file
%   addCustomTFSFSource('YeeCells', [0 0 0 100 100 100], 'SpaceTimeFile', ...
%       'data', 'Symmetries', [1 1 0]) will tell Trogdor to load E/H fields from
%       the file 'data' with a total-field scattered-field box.  The hint
%       'Symmetries' indicates that the source wave is constant along the x and
%       y directions; this may result in a substantial savings of memory and
%       runtime.
%
%   A custom TFSF source is implemented using a data request file.  Trogdor will
%   run once to generate an m-file containing an ordered list of points and 
%   timesteps at which to provide E and H fields; information about the material
%   in each cell of the grid around the TFSF boundary will be provided as well.
%   After the E and H fields have been written into the data file, Trogdor will
%   load them at the TFSF boundary and run the simulation.  This technique is
%   useful for analytical field propagation.
%
%   Usage: addCustomTFSFSource(named parameters)
%
%   Named parameters:
%       YeeCells    The cells making up the total-field region of the grid;
%                   [x0 y0 z0 x1 y1 z1] specifies all cells (x, y, z)
%                   where x0 <= x <= x1, y0 <= y <= y1, z0 <= z <= z1.
%                   (required)
%       Duration    The range of timesteps on which to source fields; [t0 t1]
%                   will source on timesteps t such that t0 <= t <= t1.  Using
%                   multiple rows specifies multiple ranges of timesteps.
%                   (default: [0, numT - 1] for a simulation with numT timesteps)
%       SpaceTimeFile   Name of a file that will contain the fields for each
%                   cell at each timestep.  Trogdor's data request file will
%                   be named [SpaceTimeFile, '.m'].  The fields in the file must 
%                   use SI units (e in V/m, h in A/m).  (required)
%       Symmetries  A three-vector indicating the symmetries of the source wave.
%                   Omitting symmetries may increase simulation time and
%                   require more data in the SpaceTimeFile.  Examples:
%                   [1 0 0]     source wave has same value for all x
%                   [0 1 1]     source wave has same value for all y and z
%       OmitSide    An axis-aligned unit vector or cell array of axis-aligned
%                   unit vectors.  The unit vector [1 0 0] specifies the +X
%                   side of the TFSF boundary and instructs Trogdor to not
%                   add or subtract electromagnetic fields on that side.
%                   (default: empty cell array)
%
%   Example:
%
%   addCustomTFSFSource('YeeCells', [0 0 0 100 100 100], 'SpaceTimeFile', ...
%       'data', 'Symmetries', [0 1 0], 'OmitSide', {[0 0 1], [0 -1 0]});
%
%   Source data into the box [0 0 0 100 100 100] but omit the correction on the
%   +Z and -Y faces.

grid = t5.TrogdorSimulation.instance().currentGrid();

% note that there is no X.MaskFile because it will ALWAYS be requested along
% with a TimeFile.  Likewise there is no Field because all E and H fields
% will be requested!
X.YeeCells = [];
X.Duration = [];
X.SpaceTimeFile = [];
X.OmitSide = [];
X.Symmetries = [0 0 0];
X = parseargs(X, varargin{:});

if size(X.YeeCells) ~= [1 6]
    error('YeeCells must have six columns.');
end

if isempty(X.Duration)
    X.Duration = [0, t5.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end

if isempty(X.SpaceTimeFile)
    error('SpaceTimeFile not provided.');
end

if any(X.Symmetries ~= 0 & X.Symmetries ~= 1)
    error('Symmetries must be a 3-vector of 1s or 0s.');
end

obj = struct;
obj.type = 'TFSFSource';
obj.yeeCells = X.YeeCells;
obj.duration = X.Duration;
obj.spaceTimeFile = X.SpaceTimeFile;
obj.symmetries = X.Symmetries;
obj.omitSides = {};

% Validate omitted sides, kind of
if ~isempty(X.OmitSide)
    if ~iscell(X.OmitSide)
        obj.omitSides = {X.OmitSide};
    else
        obj.omitSides = X.OmitSide;
    end
end

grid.CustomTFSFSources = {grid.CustomTFSFSources{:}, obj};
