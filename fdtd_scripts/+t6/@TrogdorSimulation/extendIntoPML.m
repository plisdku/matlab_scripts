function v = extendIntoPML(obj, vertices)

if size(vertices, 2) ~= 3
    error('Vertex array must be Nx3');
end

v = vertices;

outerBounds = obj.OuterBounds;
innerBounds = obj.NonPMLBounds;

for xyz = 1:3
if innerBounds(xyz) ~= innerBounds(xyz+3)
    iFloor = vertices(:,xyz) <= innerBounds(xyz);
    iCeil = vertices(:,xyz) >= innerBounds(xyz+3);
    
    v(iFloor,xyz) = outerBounds(xyz);
    v(iCeil,xyz) = outerBounds(xyz+3);
end
end





