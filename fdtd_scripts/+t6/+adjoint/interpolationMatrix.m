function matrix = interpolationMatrix(xIn, xOut)
% interpolationMatrix  Create a sparse linear interpolation matrix
%
% m = interpolationMatrix(xIn, xOut) returns a matrix of size [nOut nIn]
% where nOut = length(xOut) and nIn = length(xIn).  xOut and xIn are the
% positions of samples in an output vector and input vector.

numIn = length(xIn);
numOut = length(xOut);

if numIn == 1 && numOut == 1
    matrix = 1;
    
    if ~isequal(xIn, xOut)
        %warning('Faking interpolation matrix, please check this');
    end
    
    return;
end

% Clamp the output positions to not exceed the provided input
% positions.  Extrapolation to outside values will consequently be by
% clamping of course.  I don't need anything else right now.
xOut(xOut < xIn(1)) = xIn(1);
xOut(xOut > xIn(end)) = xIn(end);

dx = (xIn(end)-xIn(1))/(length(xIn)-1);

leftCell = 1+floor( (xOut-xIn(1))/dx ); % these are all valid indices
rightCell = min(leftCell+1, numIn); % clamp to valid indices

distLeft = xOut - xIn(leftCell);
distRight = dx - distLeft;

weightLeft = distRight/dx;
weightRight = 1 - weightLeft;

iLeft = 1:numOut;
jLeft = leftCell;

iRight = 1:numOut;
jRight = rightCell;

matrix = sparse(iLeft, jLeft, weightLeft, numOut, numIn) + ...
    sparse(iRight, jRight, weightRight, numOut, numIn);


