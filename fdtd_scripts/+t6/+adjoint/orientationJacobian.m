function Do = orientationJacobian(node, param, orientation)
% Do = orientationJacobian(node, param, orientation)
%
% node is a t6.model.Node
%
% Indexing: x y z i j fieldXYZ parameter
%

%node = group;
%param = 0;

dodv = orientationSensitivity(); % size is [numVerts, i3, j3, octants] -> [dx dy dz]

%orientation = o; % size: [x y z i j eps_i eps_j]
Do = zeros([size(orientation), length(param)]);
jacobian = node.jacobian(param);
numVerts = size(dodv,1);
numParams = size(param);

octantE = [2 3 5];

for vv = 1:numVerts
for ii = 1:3;
for jj = 1:3;
for fieldXYZ = 1:3
    octant = octantE(fieldXYZ);
if isstruct(dodv{vv, ii, jj, octant})
    
    ddv = dodv{vv, ii, jj, octant};
    
    for run = 1:size(ddv.indices,1)
        cells = ddv.indices(run,:) + 1;
        
        for pp = 1:numParams
            row = 3*(vv-1)+1;
            dotProd = dot(ddv.values(run,:), jacobian(row:row+2,pp));
            
            if ndims(Do) == 7
                Do(cells(1):cells(4), cells(2):cells(5), cells(3):cells(6),...
                    ii, jj, fieldXYZ, fieldXYZ, pp) = ...
                    Do(cells(1):cells(4),cells(2):cells(5),cells(3):cells(6),...
                        ii, jj, fieldXYZ, fieldXYZ, pp) + dotProd;
            elseif ndims(Do) == 6
                Do(cells(1):cells(4), cells(2):cells(5), cells(3):cells(6),...
                    ii, jj, fieldXYZ, pp) = ...
                    Do(cells(1):cells(4),cells(2):cells(5),cells(3):cells(6),...
                        ii, jj, fieldXYZ, pp) + dotProd;
            end
        end
    end
    
end
end
end
end
end


