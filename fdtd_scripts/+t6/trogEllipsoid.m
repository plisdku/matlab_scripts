function trogEllipsoid(rect, varargin);
% trogEllipsoid       Plot 3D outline of an ellipsoid
%   trogEllipsoid(rect) where rect bounds an ellipsoid in trogdor format, i.e.
%       [x1 y1 z1 x2 y2 z2],
%   will draw a wireframe ellipsoid inscribed in the rect with corners at
%       [x1 y1 z1]   and    [x2+1 y2+1 z2+1]
%   (since trogdor rectangle coordinates are measured in discrete indices).
%
%   trogEllipsoid(rect, 'Color', 'k') or
%   trogEllipsoid(rect, 'Color', [0 0 0]) will plot a black rectangle; other
%   Matlab colors are accepted.  Likewise any argument accepted by the line
%   command will be accepted for trogRect.
%
%   See also: trogRect, line
%
%   version 4.5
%   July 29, 2008

colors = [0 0 1];

if (nargin > 1)
    colors = varargin{1};
end

rect(4:6) = rect(4:6)+1;

radius = 0.5 * (rect(4:6)-rect(1:3));
center = 0.5 * (rect(1:3) + rect(4:6));

[x, y, z] = ellipsoid(center(1), center(2), center(3), radius(1), radius(2), ...
    radius(3));

held = ishold;

if (~ishold)
    hold on;
end

mesh(x,y,z, 'FaceAlpha', 0, 'EdgeColor', colors);

if (~ishold)
    hold off;
end