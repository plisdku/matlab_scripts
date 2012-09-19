function trogRect(rect, varargin)
% trogRect       Plot a 3D bounding rect
%   trogRect(rect) where rect is in trogdor format, i.e.
%       [x1 y1 z1 x2 y2 z2],
%   will draw the edges of a 3D rectangle with corners at
%       [x1 y1 z1]   and    [x2+1 y2+1 z2+1]
%   (since trogdor rectangle coordinates are measured in discrete indices).
%
%   trogRect('Color', 'k') or
%   trogRect('Color', [0 0 0]) will plot a black rectangle; other Matlab
%   colors are accepted.  Likewise any argument accepted by the line
%   command will be accepted for trogRect.
%
%   See also: trogEllipsoid, getOutputBounds, line
%
%   version 4.5
%   July 29, 2008

colors = [0 0 1];

if (nargin > 1)
    colors = varargin{1};
end

extras = {};
if (nargin > 2)
    [extras{1:nargin-2}] = deal(varargin{2:end});
end

rect = reshape(rect, [], 6);

for row = 1:size(rect, 1)
    r = rect(row,:);
    
    r(4:6) = r(4:6);
    r = r + [0.1 0.1 0.1 -0.1 -0.1 -0.1];

    line(r([1,4]), r([2,2]), r([3,3]), 'Color', colors, extras{:});
    line(r([1,4]), r([2,2]), r([6,6]), 'Color', colors, extras{:});
    line(r([1,4]), r([5,5]), r([3,3]), 'Color', colors, extras{:});
    line(r([1,4]), r([5,5]), r([6,6]), 'Color', colors, extras{:});

    line(r([1,1]), r([2,5]), r([3,3]), 'Color', colors, extras{:});
    line(r([4,4]), r([2,5]), r([3,3]), 'Color', colors, extras{:});
    line(r([1,1]), r([2,5]), r([6,6]), 'Color', colors, extras{:});
    line(r([4,4]), r([2,5]), r([6,6]), 'Color', colors, extras{:});

    line(r([1,1]), r([2,2]), r([3,6]), 'Color', colors, extras{:});
    line(r([4,4]), r([2,2]), r([3,6]), 'Color', colors, extras{:});
    line(r([1,1]), r([5,5]), r([3,6]), 'Color', colors, extras{:});
    line(r([4,4]), r([5,5]), r([3,6]), 'Color', colors, extras{:});
end