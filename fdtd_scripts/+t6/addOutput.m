function addOutput(filename, fields, varargin)
%addOutput Add an electromagnetic field or current output to the simulation
%   addOutput('electricFields', 'ex ey ez', 'YeeCells', [0 0 0 100 100 100])
%       will instruct Trogdor to save the ex, ey and ez fields in a binary file
%       named 'electricFields' on every timestep.
%   
%   Usage: addOutput(filename, fields, named parameters)
%   
%   Valid fields:
%       'ex', 'ey', 'ez'        electric fields
%       'hx', 'hy', 'hz'        magnetic fields
%       'jx', 'jy', 'jz'        electric currents
%       'kx', 'ky', 'kz'        magnetic currents
%   E and H fields may be saved in the same file, and J and K fields may be
%   saved in the same file, but at present currents and fields must be saved in
%   separate files.
%
%   Named parameters:
%       YeeCells        The region of the grid to save; [x0 y0 z0 x1 y1 z1] will
%                       save all cells (x, y, z) satisfying x0 <= x <= x1,
%                       y0 <= y <= y1, z0 <= z <= z1.  Multiple-row arrays will
%                       cause Trogdor to save multiple regions in one file.
%                       (required)
%       Duration        The range of timesteps on which to save data; [t0 t1]
%                       will save all timesteps t such that t0 <= t <= t1.
%                       Multiple-row arrays will cause Trogdor to save multiple
%                       ranges of timesteps; to save only timesteps 100 and 144,
%                       use addOutput('Duration', [100 100; 144 144]).
%                       (default: all timesteps)
%       Stride          Spatial sampling period.  Set to [2 2 2] to save every
%                       second cell in X, Y and Z (cutting file size by 8).  If
%                       there are multiple rows in YeeCells, the same Stride
%                       will be used for every region; if Stride has the same
%                       number of rows as YeeCells, then each region will have
%                       its own spatial sampling period.  (default: [1 1 1])
%       Period          Temporal sampling period.  Set to 10 to save every tenth
%                       timestep of each range in Duration.  If there are 
%                       multiple rows in Duration, then each range of timesteps
%                       will have the same Period; if Period has the same
%                       number of rows as Duration, then each Duration will
%                       have its own Period.  (default: 1)
%       InterpolationPoint
%                       A point between [0 0 0] and [1 1 1] at which to
%                       sample E and H fields.  The electromagnetic fields in
%                       FDTD are calculated at different spatial and temporal
%                       positions; using an InterpolationPoint will simulate
%                       measuring the fields all at one location in space by
%                       trilinear interpolation.  This will NOT resample in
%                       time.  Applicable only to E and H fields.
%                       Use this to measure E and H at the same point!
%                       (default: unused)
%           
grid = t6.TrogdorSimulation.instance().currentGrid();

X.YeeCells = []; % [x y z x y z]
X.Duration = [];  % [first last]
X.Stride = []; % scalar
X.Period = []; % scalar
X.InterpolationPoint = []; % [x y z] from 0 to 1

X = parseargs(X, varargin{:});

if length(X.YeeCells) == 0
    error('YeeCells required.');
end

if size(X.YeeCells, 2) ~= 6
    error('YeeCells must have six columns.');
end
if length(X.Duration) == 0
    X.Duration = [0, t6.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end
if length(X.Stride) == 0
    X.Stride = ones(size(X.YeeCells, 1), 3);
elseif size(X.Stride, 2) ~= 3
    error('Stride must have three columns (x, y, z).');
elseif size(X.Stride, 1) == 1
    X.Stride = repmat(X.Stride, size(X.YeeCells, 1), 1);
elseif size(X.Stride, 1) ~= size(X.YeeCells, 1)
    error('Stride must have as many rows as YeeCells.');
end
if length(X.Period) == 0
    X.Period = ones(size(X.Duration, 1), 1);
elseif size(X.Period, 1) == 1
    X.Period = repmat(X.Period, size(X.Duration, 1), 1);
elseif size(X.Period, 1) ~= size(X.Duration, 1)
    error('Period must have as many rows as Duration.');
end

if length(X.InterpolationPoint) ~= 0
    if length(X.InterpolationPoint) ~= 3
        error('Please provide a 3D interpolation point.');
    end
    if any(X.InterpolationPoint < 0) || any(X.InterpolationPoint > 1)
        error('Interpolation point must be between [0 0 0] and [1 1 1].');
    end
end

obj = struct;
obj.type = 'Output';
obj.fields = fields;
obj.filename = filename;
obj.yeeCells = X.YeeCells; % validated
obj.duration = X.Duration; % validated
obj.stride = X.Stride; % validated
obj.period = X.Period; % validated

if length(X.InterpolationPoint) ~= 0
    obj.interpolationPoint = X.InterpolationPoint;
end

grid.Outputs = {grid.Outputs{:}, obj};