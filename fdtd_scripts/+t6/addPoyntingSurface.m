function addPoyntingSurface(fileName, varargin)
%addPoyntingSurface Add an output that saves E and H on the faces of a box
%   addPoyntingSurface('poynting', 'Bounds', [0 0 0 20 20 20]) will
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
%       Bounds          The region of the simulation to save, in real units;
%                       [x0 y0 z0 x1 y1 z1] will determine the Yee cell rect
%                       [m0 n0 p0 m1 n1 p1] to save, suitably for the grid
%                       resolution. (required)
%       Duration        The range of timesteps on which to save data; [t0 t1]
%                       will save all timesteps t such that t0 <= t <= t1.
%                       Multiple-row arrays will cause Trogdor to save multiple
%                       ranges of timesteps; to save only timesteps 100 and 144,
%                       use addOutput('Duration', [100 100; 144 144]).
%                       (default: all timesteps)
%       CutoffFrequency The highest angular frequency expected to be
%                       measured.  The sampling period of the output file
%                       will be the largest possible integer number of
%                       timesteps that still resolves the cutoff frequency.
%                       Thus the sampling rate will be at least twice the
%                       cutoff frequency according to the Nyquist sampling
%                       theorem.  CutoffFrequency is an alternative to Period.
%                       (default: unused)
%       Period          Temporal sampling period.  Set to 10 to save every tenth
%                       timestep of each range in Duration.  If there are 
%                       multiple rows in Duration, then each range of timesteps
%                       will have the same Period; if Period has the same
%                       number of rows as Duration, then each Duration will
%                       have its own Period.  (default: 1)
%       Sides           A row vector of six elements.  A nonzero value
%                       at position n signifies that the output should
%                       include face n.  The order of faces is low X, high
%                       X, low Y, high Y, low Z, high Z.
%                       (default: [1 1 1 1 1 1])

X.Bounds = [];
X.Duration = [];
X.Period = [];
X.CutoffFrequency = [];
X.Sides = [1 1 1 1 1 1];
X = parseargs(X, varargin{:});

t6.addSurfaceOutput(fileName, 'ex ey ez hx hy hz', ...
    'Bounds', X.Bounds, 'Duration', X.Duration, ...
    'Period', X.Period, 'CutoffFrequency', X.CutoffFrequency, ...
    'Sides', X.Sides);
