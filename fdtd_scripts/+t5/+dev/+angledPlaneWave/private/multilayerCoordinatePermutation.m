function [ijk, ijkInv] = multilayerCoordinatePermutation(afp, direction)

ijk = [1 2 3];
ijkInv = [1 2 3];
if afp.sampleHalfCells(1) ~= afp.sampleHalfCells(4)
    ijk = [2 3 1];
    ijkInv = [3 1 2];
elseif afp.sampleHalfCells(2) ~= afp.sampleHalfCells(5)
    ijk = [3 1 2];
    ijkInv = [2 3 1];
elseif afp.sampleHalfCells(3) ~= afp.sampleHalfCells(6)
    ijk = [1 2 3];
    ijkInv = [1 2 3];
else
    warning('Only one material present.  Unknown orientation.');
    
    [nothing, dir] = max(abs(direction));
    if dir == 1 % x
        ijk = [2 3 1];
        ijkInv = [3 1 2];
    elseif dir == 2 % y
        ijk = [3 1 2];
        ijkInv = [2 3 1];
    else
        ijk = [1 2 3];
        ijkInv = [1 2 3];
    end
end
