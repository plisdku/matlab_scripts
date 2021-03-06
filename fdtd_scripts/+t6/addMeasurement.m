function addMeasurement(varargin)
%addMeasurement Add a measurement to optimize
%   addMeasurement('Fields', 'ex ey ez', 'Bounds', [0 0 0 100 100 100], ...
%       'Filters', f)
%       will prepare Trogdor to calculate the function f and its derivative
%       df/dFields from a forward FDTD simulation, and to generate a current
%       source for adjoint simulations.
%   
%   Usage: addMeasurement(named parameters)
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
%       Fields          Which electromagnetic fields to save
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
%       Filters         The filters for the merit function
%       Kernel          The kernel for the merit function (optional)
%       Timesteps       The range of timesteps on which to save data; [t0 t1]
%                       will save all timesteps t such that t0 <= t <= t1.
%                       Multiple-row arrays will cause Trogdor to save multiple
%                       ranges of timesteps; to save only timesteps 100 and 144,
%                       use addOutput('Timesteps', [100 100; 144 144]).
%                       Timesteps range from 0 to numTimesteps-1.
%                       (default: all timesteps)
%       Stride          Spatial sampling period.  Set to [2 2 2] to save every
%                       second cell in X, Y and Z (cutting file size by 8).  If
%                       there are multiple rows in YeeCells, the same Stride
%                       will be used for every region; if Stride has the same
%                       number of rows as YeeCells, then each region will have
%                       its own spatial sampling period.  (default: [1 1 1])
%       Period          Temporal sampling period.  Set to 10 to save every tenth
%                       timestep of each range in Timesteps.  If there are 
%                       multiple rows in Timesteps, then each range of timesteps
%                       will have the same Period; if Period has the same
%                       number of rows as Timesteps, then each Timesteps will
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

sim = simulation();

X.Fields = '';
X.YeeCells = []; % [x y z x y z]
X.Bounds = [];
X.Function = [];
X.Filters = [];
X.Kernel = [];
X.Timesteps = [];  % [first last]
X.Stride = []; % scalar
X.Period = []; % scalar
X.CutoffFrequency = [];
X.InterpolationPoint = []; % [x y z] from 0 to 1
X.Mode = '';

X = parseargs(X, varargin{:});

validateYeeCellsAndBounds(X);

fieldTokens = mySortTokens(tokenizeFields(X.Fields, 'd e b h j m'));

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
if ~isempty(X.Bounds)
    X.YeeCells = sim.boundsToYee(X.Bounds, fieldTokens);
end

if isempty(X.Timesteps)
    X.Timesteps = [0, sim.NumT-1];
elseif size(X.Timesteps, 2) ~= 2
    error('Timesteps must have two columns (first and last timestep).');
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
        X.Period = floor(0.5*2*pi./X.CutoffFrequency/sim.Dt);
        X.Period(X.Period < 1) = 1;
    else
        X.Period = ones(size(X.Timesteps, 1), 1);
    end
end

if size(X.Period, 1) == 1
    X.Period = repmat(X.Period, size(X.Timesteps, 1), 1);
elseif size(X.Period, 1) ~= size(X.Timesteps, 1)
    error('Period must have as many rows as Timesteps.');
end

if ~isempty(X.InterpolationPoint)
    if length(X.InterpolationPoint) ~= 3
        error('Please provide a 3D interpolation point.');
    end
    if any(X.InterpolationPoint < 0) || any(X.InterpolationPoint > 1)
        error('Interpolation point must be between [0 0 0] and [1 1 1].');
    end
end

if ~isempty(X.Mode)
    if ~strcmpi(X.Mode, 'forward') && ~strcmpi(X.Mode, 'Adjoint')
        error('Mode must be forward or adjoint if specified');
    end
end

if ~isempty(X.Function)
    X.Filters = X.Function;
    warning('Function is deprecated and should be replaced with Filters');
end

obj = struct;
obj.type = 'Output';
obj.fields = X.Fields;
obj.filters = X.Filters;
obj.kernel = X.Kernel;
obj.filename = sprintf('measurement%i', length(sim.Grid.Measurements)+1);
obj.yeeCells = X.YeeCells; % validated
obj.bounds = X.Bounds; % validated
obj.timesteps = X.Timesteps; % validated
obj.stride = X.Stride; % validated
obj.period = X.Period; % validated
obj.mode = X.Mode;

if ~isempty(X.InterpolationPoint)
    obj.interpolationPoint = X.InterpolationPoint;
end

%if ~isempty(sim.Grid.Measurement)
%    warning('Overwriting measurement');
%end

sim.Grid.Measurements{end+1} = obj;

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


