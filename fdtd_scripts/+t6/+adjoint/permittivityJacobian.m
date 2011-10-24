function [num, den] = permittivityJacobian(node, param, epsNumer, ...
    epsDenom, varargin)
% [num, den] = permittivityJacobian(node, param, epsNumer, epsDenom)
%
% [num, den] = permittivityJacobian(node, param, epsNumer, epsDenom, 
%   vertexNumbers)
%
% [num, den] = permittivityJacobian(node, param, epsNumer, epsDenom, 
%   vertexNumbers, freeDirections)
%
% node is a t6.model.Node
%
% Indexing: [x y z fieldXYZ lag paramIndex]
%

%node = group;
%param = parameters(pp);
%epsNumer = num;
%epsDenom = den;

[coeffs, verts] = adjoint.readDeps();

if nargin > 4
    verts = varargin{1};
end

freeDirections = 1:size(coeffs, 2);
if nargin > 5
    freeDirections = intersect(varargin{2}, 1:size(coeffs,2));
end

jacobian = node.jacobian(param);

num = zeros([size(epsNumer), length(param)]);
den = zeros([size(epsDenom), length(param)]);

for vv = verts
for freeDirection = freeDirections % 1:3 at most (xyz)
    
    dDBdCoord = zeros(size(epsDenom)); % sensitivity w.r.t. (e.g.) v3y
    
    if isstruct(coeffs{vv,freeDirection})
    for fieldXYZ = 1:length(coeffs{vv,freeDirection}.tensor)
    if isstruct(coeffs{vv,freeDirection}.tensor{fieldXYZ, fieldXYZ})
        Ddb = coeffs{vv, freeDirection}.tensor{fieldXYZ, fieldXYZ}.DB;

        for lag = 1:length(Ddb)
        for pp = 1:size(Ddb{lag}.positions,1)
            cellX = Ddb{lag}.positions(pp,1) + 1;
            cellY = Ddb{lag}.positions(pp,2) + 1;
            %cellZ = Ddb{lag}.positions(pp,3) + 1;
            cellZ = 1;

            dDBdCoord(cellX, cellY, cellZ, fieldXYZ, fieldXYZ, lag) = ...
                dDBdCoord(cellX, cellY, cellZ, fieldXYZ, fieldXYZ, lag) + ...
                Ddb{lag}.coefficients(pp);
        end
        end
    end
    end
    end
    % Jacobian: num = dDBdCoord * jacobian
    for pp = 1:length(param)
    %for jj = 1:size(jacobian,1);
        row = 3*(vv-1) + freeDirection;
        den(:,:,:,:,:,:,pp) = den(:,:,:,:,:,:,pp) ...
            + dDBdCoord * full(jacobian(row,1));
    end
    %end
end
end


for vv = verts % v1x v1y v1z v2x v2y v2z ...
for freeDirection = 1:size(coeffs,2) % 1:3 at most (xyz)
    
    dEHdCoord = zeros(size(epsNumer)); % sensitivity w.r.t. (e.g.) v3y
    
    if isstruct(coeffs{vv,freeDirection})
    for fieldXYZ = 1:length(coeffs{vv,freeDirection}.tensor)
    if isstruct(coeffs{vv,freeDirection}.tensor{fieldXYZ, fieldXYZ})
        Deh = coeffs{vv, freeDirection}.tensor{fieldXYZ, fieldXYZ}.EH;

        for lag = 1:length(Deh)
        for pp = 1:size(Deh{lag}.positions,1)
            cellX = Deh{lag}.positions(pp,1) + 1;
            cellY = Deh{lag}.positions(pp,2) + 1;
            %cellZ = Ddb{lag}.positions(pp,3) + 1;
            cellZ = 1;

            dEHdCoord(cellX, cellY, cellZ, fieldXYZ, fieldXYZ, lag) = ...
                dEHdCoord(cellX, cellY, cellZ, fieldXYZ, fieldXYZ, lag) + ...
                Deh{lag}.coefficients(pp);
        end
        end
    end
    end
    end
    
    % Jacobian: num = dDBdCoord * jacobian
    for pp = 1:length(param)
    %for jj = 1:size(jacobian,1);
        row = 3*(vv-1) + freeDirection;
        num(:,:,:,:,:,:,pp) = num(:,:,:,:,:,:,pp) ...
            + dEHdCoord * full(jacobian(row,1));
    end
    %end
end
end
