function axisString = makeAxisString(axisVector)
%makeAxisString Return a string representation of an axis-oriented unit vector
%   This function is internal to Trogdor.

if length(axisVector) ~= 3
    error('Axis vector must be a three-dimensional unit vector like [1 0 0].');
end

if axisVector(1) > 0
    axisString = 'x';
elseif axisVector(1) < 0
    axisString = '-x';
elseif axisVector(2) > 0
    axisString = 'y';
elseif axisVector(2) < 0
    axisString = '-y';
elseif axisVector(3) > 0
    axisString = 'z';
else
    axisString = '-z';
end

