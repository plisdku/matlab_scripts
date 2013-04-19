function addOutput(filename, fields, varargin)
%addOutput Add an electromagnetic field or current output to the simulation
%   addOutput('electricFields', 'ex ey ez', 'Bounds', [0 0 0 100 100 100])
%       will instruct Trogdor to save the ex, ey and ez fields in a file
%       named 'electricFields' on every timestep.
%   
%   Usage: addOutput(filename, fields, named parameters)
%          addOutput(a, b, c)
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
%       Duration        The range of time over which data will be saved.
%                       [t0 t1] will save enough timesteps to find values
%                       of fields in the time interval [t0 t1].
%       Timesteps       The range of timesteps on which to save data; [n0 n1]
%                       will save all timesteps n such that n0 <= n <= n1.
%                       Multiple-row arrays will cause Trogdor to save multiple
%                       ranges of timesteps; to save only timesteps 100 and 144,
%                       use addOutput('Timesteps', [100 100; 144 144]).
%                       Timesteps range from 0 to numTimesteps-1.
%                       (default: all timesteps)
%       Frequency       Frequency or frequencies of on-the-fly Fourier
%                       transform (default: unused)
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

X.YeeCells = []; % [x y z x y z]
X.Bounds = [];
X.Timesteps = [];  % [first last]
X.Frequency = [];
X.Duration = []; % [first last]
X.Period = []; % scalar
X.CutoffFrequency = [];
X.InterpolationPoint = []; % [x y z] from 0 to 1
X.Mode = '';

X = parseargs(X, varargin{:});

validateYeeCellsAndBounds(X);

fieldTokens = mySortTokens(tokenizeFields(fields, 'd e b h j m'));

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
if ~isempty(X.Bounds)
    X.YeeCells = sim.boundsToYee(X.Bounds, fieldTokens);
end

% Make sure Frequency does not coexist with Timesteps or Duration.
% It can coexist with Period and CutoffFrequency I guess, sure, since those
% can be used to save some computation time in FDTD.

if ~isempty(X.Frequency) && (~isempty(X.Timesteps) || ~isempty(X.Duration))
    error('Frequency option cannot coexist with Timesteps or Duration');
end

% If we obtained Duration and not Timesteps, set the Timesteps
% appropriately

if ~isempty(X.Timesteps) && size(X.Timesteps, 2) ~= 2
    error('Timesteps must have two columns (first and last timestep)');
end

if ~isempty(X.Duration) && size(X.Duration, 2) ~= 2
    error('Duration must have two columns (beginning and end time)');
end

if ~isempty(X.Duration)
    X.Timesteps = sim.timeToTimesteps(X.Duration, fieldTokens);
end

if isempty(X.Timesteps)
    X.Timesteps = [0, sim.NumT-1];
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


obj = struct;
obj.type = 'Output';
obj.fields = fields;
obj.filename = filename;
obj.yeeCells = X.YeeCells; % validated
obj.bounds = X.Bounds; % validated
obj.timesteps = X.Timesteps; % validated
obj.frequency = X.Frequency;
obj.duration = X.Duration; % validated
obj.period = X.Period; % validated
obj.mode = X.Mode;

if ~isempty(X.InterpolationPoint)
    obj.interpolationPoint = X.InterpolationPoint;
end

sim.Grid.Outputs{end+1} = obj;

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


