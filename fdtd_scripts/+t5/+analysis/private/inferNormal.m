
function normal = inferNormal(positions, regionNumber)

% 1.  Obtain entire bounding box
boundingBox = [inf inf inf -inf -inf -inf];

for rr = 1:size(positions, 1)
    for xyz = 1:3
        boundingBox(xyz) = min([boundingBox(xyz), positions{rr,xyz}]);
        boundingBox(xyz+3) = max([boundingBox(xyz+3), positions{rr,xyz}]);
    end
end

normal = [0 0 0]';
for xyz = 1:3
    if numel(positions{regionNumber,xyz}) == 1
        if positions{regionNumber,xyz} == boundingBox(xyz) && ...
            positions{regionNumber,xyz} ~= boundingBox(xyz+3)
            normal(xyz) = -1;
            return
        elseif positions{regionNumber,xyz} == boundingBox(xyz+3) && ...
            positions{regionNumber,xyz} ~= boundingBox(xyz)
            normal(xyz) = 1;
            return;
        end
    end
end
error('Cannot infer normal vector');

