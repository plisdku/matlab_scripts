function Df = fillFactorJacobian(node, param, fillFactors)
% Df = fillFactorJacobian(node, param, fillFactors)
%
% node is a t6.model.Node
%
% Indexing: [x y z material fieldXYZ parameter]

%node = group;
%param = 0;
%fillFactors = f; % size: [x y z material fieldXYZ]

%if ndims(fillFactors) ~= 6
%    error('Fill factor array does not seem to have all octants');
%end

dfdv = fillFactorSensitivity(); % size is [numOctants numMaterials numVertices]

octantsE = [2 3 5];

Df = zeros([size(fillFactors), length(param)]);
jacobian = node.jacobian(param);
numVerts = size(dfdv,3);
numMaterials = size(dfdv,2);
numParams = size(param);

for mm = 1:numMaterials
for vv = 1:numVerts
for fieldXYZ = 1:3
    octant = octantsE(fieldXYZ);
if isstruct(dfdv{octant,mm,vv})
    
    ddv = dfdv{octant,mm,vv};
    for run = 1:size(ddv.indices,1)
        cells = ddv.indices(run,:) + 1;
        
        for pp = 1:numParams
            row = 3*(vv-1)+1;
            dotProd = dot(ddv.values(run,:), jacobian(row:row+2,pp));
            
            if ndims(fillFactors) == 6
                Df(cells(1):cells(4), cells(2):cells(5), cells(3):cells(6),...
                    mm, fieldXYZ, fieldXYZ, pp) = ...
                    Df(cells(1):cells(4),cells(2):cells(5),cells(3):cells(6),...
                        mm, fieldXYZ, fieldXYZ, pp) + dotProd;
            elseif ndims(fillFactors) == 5
                Df(cells(1):cells(4), cells(2):cells(5), cells(3):cells(6),...
                    mm, fieldXYZ, pp) = ...
                    Df(cells(1):cells(4),cells(2):cells(5),cells(3):cells(6),...
                        mm, fieldXYZ, pp) + dotProd;
            end
        end
    end
    
end
end
end
end
