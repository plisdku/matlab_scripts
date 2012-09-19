function yeeCells = boundsToYee(obj, bounds, fieldTokens)
% yeeCells = boundsToYee(bounds, fieldTokens)
%
% bounds should be a real-valued six-vector [x0 y0 z0 x1 y1 z1]
% fieldTokens should be a cell array of field names, e.g. {'ex', 'jy', 'mz'}

import t6.*

if ~iscell(fieldTokens)
    fieldTokens = {fieldTokens};
end

offset = t6.fieldOffset(fieldTokens{1}) .* [obj.Dxyz 1];

yeeCells = obj.yeeCells(bsxfun(@minus, bounds, [offset(1:3) offset(1:3)]));
for ff = 2:numel(fieldTokens)
    offset = t6.fieldOffset(fieldTokens{ff}) .* [obj.Dxyz 1];
    for rr = 1:size(bounds, 1)
        yeeCells(rr,:) = rectUnion(yeeCells(rr,:), ...
            obj.yeeCells(bounds(rr,:) - [offset(1:3) offset(1:3)]));
    end
end
