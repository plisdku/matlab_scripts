function addPoyntingOutput(fileName, varargin)
%addPoyntingOutput Add an output that saves ex, ey, ez, hx, hy, and hz
%   addPoyntingOutput('poynting', 'Bounds', [0 0 0 0 20 20]) will instruct
%   Trogdor to save the E and H fields in a binary file named 'poynting' on
%   every timestep.
%
%   This function differs from addOutput essentially in specifying the six
%   E and H field components by default.
%
%   Usage: addPoyntingOutput(filename, named parameters)
%
%   Named parameters:
%       Bounds          The region of the simulation to save, in real units;
%                       [x0 y0 z0 x1 y1 z1] will determine the Yee cell rect
%                       [m0 n0 p0 m1 n1 p1] to save, suitably for the grid
%                       resolution. (required)
%       Timesteps       The range of timesteps on which to save data; [t0 t1]
%                       will save all timesteps t such that t0 <= t <= t1.
%                       Multiple-row arrays will cause Trogdor to save multiple
%                       ranges of timesteps; to save only timesteps 100 and 144,
%                       use addOutput('Timesteps', [100 100; 144 144]).
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
%                       timestep of each range in Timesteps.  If there are 
%                       multiple rows in Timesteps, then each range of timesteps
%                       will have the same Period; if Period has the same
%                       number of rows as Timesteps, then each Timesteps will
%                       have its own Period.  (default: 1)

X.Bounds = [];
X.Timesteps = [];
X.Period = [];
X.CutoffFrequency = [];
X = parseargs(X, varargin{:});

t6.addOutput(fileName, 'ex ey ez hx hy hz', ...
    'Bounds', X.Bounds, 'Timesteps', X.Timesteps, ...
    'Period', X.Period, 'CutoffFrequency', X.CutoffFrequency);


