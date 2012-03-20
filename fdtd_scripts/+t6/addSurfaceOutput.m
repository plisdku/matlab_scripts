function addSurfaceOutput(fileName, fields, varargin)
%addSurfaceOutput Add an output that saves fields on the faces of a box.
%   addSurfaceOutput('surf', 'ex', 'Bounds', [0 0 0 20 20 20]) will
%   instruct Trogdor to save the ex field along the six faces of
%   [0 0 0 20 20 20] in a binary file named 'surf' on every timestep.
%   The six rectangles in order will be
%       [0 0 0 0 20 20]         (low X)
%       [20 0 0 20 20 20]       (high X)
%       [0 0 0 20 0 20]         (low Y)
%       [0 20 0 20 20 20]       (high Y)
%       [0 0 0 20 20 0]         (low Z)
%       [0 0 20 20 20 20]       (high Z)
%
%   Usage: addSurfaceOutput(filename, fields, named parameters)
%
%   Named parameters:
%       Bounds          The region of the simulation to save, in real units;
%                       [x0 y0 z0 x1 y1 z1] will determine the Yee cell rect
%                       [m0 n0 p0 m1 n1 p1] to save, suitably for the grid
%                       resolution.
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
X.YeeCells = [];
X.Duration = [];
X.Stride = [1 1 1];
X.Period = [];
X.CutoffFrequency = [];
X.InterpolationPoint = [];
X.Sides = [1 1 1 1 1 1];
X = parseargs(X, varargin{:});

t6.validateYeeCellsAndBounds(X);

fieldTokens = mySortTokens(t6.tokenizeFields(fields, 'd e b h'));

% If we obtained Bounds and not YeeCells, set the YeeCells appropriately
%if ~isempty(X.Bounds)
%    X.YeeCells = boundsToYee(X.Bounds, fieldTokens);
%end

bounds = [];
for side = 1:6
    xyz = 1+floor( (side-1)/2 );
    if X.Sides(side) && X.Bounds(xyz) < X.Bounds(xyz + 3)
        bounds(end+1,:) = sideOfRect(X.Bounds, side);
    end
end

if isempty(bounds)
    warning('Poynting surface omits all six sides!');
end

t6.addOutput(fileName, fields, ...
    'Bounds', bounds, 'Duration', X.Duration, ...
    'Period', X.Period, 'CutoffFrequency', X.CutoffFrequency);


function rect = sideOfRect(inRect, whichSide)

switch whichSide
    case 1 % low X
        rect = inRect([1 2 3 1 5 6]);
    case 2 % high X
        rect = inRect([4 2 3 4 5 6]);
    case 3 % low Y
        rect = inRect([1 2 3 4 2 6]);
    case 4 % high Y
        rect = inRect([1 5 3 4 5 6]);
    case 5 % low Z
        rect = inRect([1 2 3 4 5 3]);
    case 6 % high Z
        rect = inRect([1 2 6 4 5 6]);
    otherwise
        error('Side of rect must range from 1 to 6');
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
