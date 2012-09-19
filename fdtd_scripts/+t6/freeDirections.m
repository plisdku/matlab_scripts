function freeDirections(varargin)

import t6.*
sim = simulation();
grid = sim.CurrentGrid;

X.YeeBounds = [0 0 0 0 0 0];
X.Directions = [0 0 0];
X = parseargs(X, varargin{:});

if ~validateRect(X.YeeBounds)
    error('Invalid rectangle.');
end

for nn = 1:length(grid.Assembly)
if strcmp(grid.Assembly{nn}.type, 'Mesh')
    for vv = 1:length(grid.Assembly{nn}.vertices)
        if all(grid.Assembly{nn}.vertices(vv,:) >= X.YeeBounds(1:3)) &&...
            all(grid.Assembly{nn}.vertices(vv,:) <= X.YeeBounds(4:6))
            grid.Assembly{nn}.vertexFreeDirections(vv,:) = X.Directions;
        end
    end
end
end