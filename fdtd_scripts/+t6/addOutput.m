function addOutput(filename, fields, varargin)
%addOutput Add an electromagnetic field or current output to the simulation
%   addOutput('electricFields', 'ex ey ez', 'Bounds', [0 0 0 100 100 100])
%       will instruct Trogdor to save the ex, ey and ez fields in a file
%       named 'electricFields' on every timestep.
%   
%   Usage: addOutput(filename, fields, named parameters)
%   
%   Valid fields:
%       'ex', 'ey', 'ez'        electric fields
%       'exx', 'exy', etc.      off-diagonal electric field subcomponents
%       'hx', 'hy', 'hz'        magnetic fields
%       'hxx', 'hxy', etc.      off-diagonal magnetic field subcomponents
%       'dx', 'dy', 'dz'        D field
%       'bx', 'by', 'bz'        B field 
%   The order of fields in the file will always be D, E, B, H.
%
%   Named parameters:
%       YeeCells        The region of the grid to save; [x0 y0 z0 x1 y1 z1] will
%                       save all cells (x, y, z) satisfying x0 <= x <= x1,
%                       y0 <= y <= y1, z0 <= z <= z1.  Multiple-row arrays will
%                       cause Trogdor to save multiple regions in one file.
%                       (YeeCells or Bounds required)
%       Bounds          The region of the simulation to save, in real units;
%                       [x0 y0 z0 x1 y1 z1] will determine the Yee cell rect
%                       [m0 n0 p0 m1 n1 p1] to save, suitably for the grid
%                       resolution.
%                       (YeeCells or Bounds required)
%       Duration        The range of timesteps on which to save data; [t0 t1]
%                       will save all timesteps t such that t0 <= t <= t1.
%                       Multiple-row arrays will cause Trogdor to save multiple
%                       ranges of timesteps; to save only timesteps 100 and 144,
%                       use addOutput('Duration', [100 100; 144 144]).
%                       Timesteps range from 0 to numTimesteps-1.
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
%       CutoffFrequency The highest angular frequency expected to be
%                       measured.  The sampling period of the output file
%                       will be the largest possible integer number of
%                       timesteps that still resolves the cutoff frequency.
%                       Thus the sampling rate will be at least twice the
%                       cutoff frequency according to the Nyquist sampling
%                       theorem.  CutoffFrequency is an alternative to Period.
%                       (default: unused)
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

import t6.*

grid = t6.TrogdorSimulation.instance().currentGrid();

X.YeeCells = []; % [x y z x y z]
X.Bounds = [];
X.Duration = [];  % [first last]
X.Stride = []; % scalar
X.Period = []; % scalar
X.CutoffFrequency = [];
X.InterpolationPoint = []; % [x y z] from 0 to 1

X = parseargs(X, varargin{:});

t6.validateYeeCellsAndBounds(X);

fieldTokens = mySortTokens(tokenizeFields(fields, 'd e b h j m'));

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
if ~isempty(X.Bounds)
    X.YeeCells = boundsToYee(X.Bounds, fieldTokens);
end

if isempty(X.Duration)
    X.Duration = [0, t6.TrogdorSimulation.instance().NumT-1];
elseif size(X.Duration, 2) ~= 2
    error('Duration must have two columns (first and last timestep).');
end
if isempty(X.Stride)
    X.Stride = ones(size(X.YeeCells, 1), 3);
elseif size(X.Stride, 2) ~= 3
    error('Stride must have three columns (x, y, z).');
elseif size(X.Stride, 1) == 1
    X.Stride = repmat(X.Stride, size(X.YeeCells, 1), 1);
elseif size(X.Stride, 1) ~= size(X.YeeCells, 1)
    error('Stride must have as many rows as YeeCells or Bounds.');
end

if ~isempty(X.Period) && ~isempty(X.CutoffFrequency)
    error('Period and CutoffFrequency cannot both be used');
end

if isempty(X.Period)
    if ~isempty(X.CutoffFrequency)
        X.Period = floor(0.5*2*pi./X.CutoffFrequency/t6.simulation().Dt);
        X.Period(X.Period < 1) = 1;
    else
        X.Period = ones(size(X.Duration, 1), 1);
    end
end

if size(X.Period, 1) == 1
    X.Period = repmat(X.Period, size(X.Duration, 1), 1);
elseif size(X.Period, 1) ~= size(X.Duration, 1)
    error('Period must have as many rows as Duration.');
end

if ~isempty(X.InterpolationPoint)
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
obj.bounds = X.Bounds; % validated
obj.duration = X.Duration; % validated
obj.stride = X.Stride; % validated
obj.period = X.Period; % validated

if ~isempty(X.InterpolationPoint)
    obj.interpolationPoint = X.InterpolationPoint;
end

grid.Outputs = {grid.Outputs{:}, obj};

end


function toks = mySortTokens(tokens)

vals = zeros(size(tokens));

for tt = 1:numel(tokens)
    switch lower(tokens{tt}(1))
        case 'd'
            vals(tt) = 1;
        case 'e'
            vals(tt) = 2;
        case 'b'
            vals(tt) = 3;
        case 'h'
            vals(tt) = 4;
        case 'j'
            vals(tt) = 5;
        case 'm'
            vals(tt) = 6;
    end
end

[unused, ii] = sort(vals);

if any(ii ~= 1:numel(ii))
    warning('Output fields have been reordered to D, E, B, H.');
end

toks = tokens(ii);        
    
end % function tokenPrecedence


