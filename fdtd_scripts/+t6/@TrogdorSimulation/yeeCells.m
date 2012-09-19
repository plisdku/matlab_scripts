function y = yeeCells(obj, bounds)
%yeeCells Calculate the rectangle of cells that encompasses the given
%bounds (real coordinates); clip to grid bounds.

import t6.*

dxyz = obj.Dxyz;
grid = obj.CurrentGrid;
origin = grid.Origin;

boundsCells = bsxfun(@times, bsxfun(@minus, bounds, [origin origin]), ...
    1./[dxyz dxyz]);

y = [floor(boundsCells(:,1:3)), ceil(boundsCells(:,4:6))];

for xyz = 1:3
    if grid.YeeCells(xyz) == grid.YeeCells(xyz+3)
        y(:, [xyz xyz+3]) = repmat(grid.YeeCells([xyz xyz+3]), size(y,1), 1);
    end
end

for rr = 1:size(y,1)
    y(rr,:) = [max(y(rr,1:3), grid.YeeCells(1:3)), ...
        min(y(rr, 4:6), grid.YeeCells(4:6))];
end


