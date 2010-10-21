function s = squish(a, dims)
% s = squish(a, dims)
% Squeeze out certain singleton dimensions

if ndims(a) > 2
    sizeOfA = size(a);
    sizeOfA(dims) = [];
    s = reshape(a, sizeOfA);
else
    s = a;
end