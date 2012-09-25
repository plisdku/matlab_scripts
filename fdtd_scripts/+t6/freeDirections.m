function freeDirections(varargin)

import t6.*
sim = simulation();
%grid = sim.CurrentGrid;

X.YeeBounds = [0 0 0 0 0 0];
X.Directions = [0 0 0];
X = parseargs(X, varargin{:});

if ~validateRect(X.YeeBounds)
    error('Invalid rectangle.');
end

for nn = 1:length(sim.Grid.Assembly)
if strcmp(sim.Grid.Assembly{nn}.type, 'Mesh')
    for vv = 1:length(sim.Grid.Assembly{nn}.vertices)
        if all(sim.Grid.Assembly{nn}.vertices(vv,:) >= X.YeeBounds(1:3)) &&...
            all(sim.Grid.Assembly{nn}.vertices(vv,:) <= X.YeeBounds(4:6))
            sim.Grid.Assembly{nn}.vertexFreeDirections(vv,:) = X.Directions;
        end
    end
end
end