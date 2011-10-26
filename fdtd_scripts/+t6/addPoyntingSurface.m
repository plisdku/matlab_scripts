function addPoyntingSurface(fileName, varargin)
%addPoyntingSurface Add an output that saves E and H on the faces of a box
%   addPoyntingSurface('poynting', 'YeeCells', [0 0 0 20 20 20]) will
%   instruct Trogdor to save the E and H fields along the six faces of
%   [0 0 0 20 20 20] in a binary file named 'poynting' on every timestep.
%   The six rectangles in order will be
%       [0 0 0 0 20 20]         (low X)
%       [20 0 0 20 20 20]       (high X)
%       [0 0 0 20 0 20]         (low Y)
%       [0 20 0 20 20 20]       (high Y)
%       [0 0 0 20 20 0]         (low Z)
%       [0 0 20 20 20 20]       (high Z)
%
%   This function differs from addSurfaceOutput essentially in specifying
%   the six E and H field components and interpolation point [0 0 0] by
%   default; it will fail if the Stride does not allow for output in each
%   corner of YeeCells (this is important for integration of the Poynting
%   vector flux through a surface).
%
%   Usage: addPoyntingSurface(filename, named parameters)
%
%   Named parameters:
%       YeeCells        The region of the grid to save; [x0 y0 z0 x1 y1 z1] will
%                       save all cells (x, y, z) satisfying x0 <= x <= x1,
%                       y0 <= y <= y1, z0 <= z <= z1.
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
%                       (default: [0 0 0])
%       Sides           A row vector of six elements.  A nonzero value
%                       at position n signifies that the output should
%                       include face n.  The order of faces is low X, high
%                       X, low Y, high Y, low Z, high Z.
%                       (default: [1 1 1 1 1 1])

X.YeeCells = [];
X.Duration = [];
X.Stride = [1 1 1];
X.Period = 1;
X.InterpolationPoint = [0 0 0];
X.Sides = [1 1 1 1 1 1];
X = parseargs(X, varargin{:});

t6.addSurfaceOutput(fileName, 'ex ey ez hx hy hz', ...
    'YeeCells', X.YeeCells, 'Duration', X.Duration, ...
    'Period', X.Period, 'InterpolationPoint', X.InterpolationPoint, ...
    'Stride', X.Stride, 'Sides', X.Sides);
