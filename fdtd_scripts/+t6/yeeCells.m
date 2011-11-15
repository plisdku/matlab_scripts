function y = yeeCells(bounds)
%yeeCells Calculate the rectangle of cells that encompasses the given
%bounds (real coordinates); clip to grid bounds.

import t6.*

dxyz = t6.sim().Dxyz;
origin = t6.grid().Origin;

b2 = bounds - [origin origin];

y = [floor(b2(1:3)./dxyz), ceil(b2(4:6)./dxyz)];

for xyz = 1:3
    if t6.grid().YeeCells(xyz) == t6.grid().YeeCells(xyz+3)
        y([xyz xyz+3]) = t6.grid().YeeCells([xyz xyz+3]);
    end
end

y = [max(y(1:3), t6.grid().YeeCells(1:3)), min(y(4:6), t6.grid().YeeCells(4:6))];

