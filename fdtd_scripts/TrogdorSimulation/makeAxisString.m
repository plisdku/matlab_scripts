% axisString    Internal function for Trogdor MatlabXML methods
%   
%   version 4.5
%   July 29, 2008
function axisString = makeAxisString(axisVector)

if axisVector(1) > 0
    axisString = 'x';
elseif axisVector(1) < 0
    axisString = '-x';
elseif axisVector(2) > 0  % Trogdor uses vertical columns so images are XY.
    axisString = '-y';
elseif axisVector(2) < 0
    axisString = 'y';
elseif axisVector(3) > 0
    axisString = 'z';
else
    axisString = '-z';
end

