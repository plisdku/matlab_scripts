function [f, dfdp] = loadFillFactors()

f = adjoint.readFillFactors('fillFactors');
Df = fillFactorSensitivity();



%%

o = readOrientation('orientations');
Do = orientationSensitivity(); % indexing: x y z i j fieldXYZ

% Assume that all vertices move in the x direction.

dodp = 0*o;

for vv = 1:size(Do, 1)
    for xyz = 1:3
        bounds = Do{vv,xyz,xyz,xyz}.bounds;
        indices = Do{vv,xyz,xyz,xyz}.indices;
        values = Do{vv,xyz,xyz,xyz}.values;
        
        if ~isempty(indices)
            ii = indices + 1;
            
            for nn = 1:size(ii,1)
            for freeDirection = 1 % CAREFUL
                dodp(ii(nn,1):ii(nn,4), ii(nn,2):ii(nn,5), ii(nn,3):ii(nn,6), ...
                    xyz, xyz, xyz) = ...
                dodp(ii(nn,1):ii(nn,4), ii(nn,2):ii(nn,5), ii(nn,3):ii(nn,6), ...
                    xyz, xyz, xyz) + values(nn,freeDirection);
            end
            end
        end
    end
end

