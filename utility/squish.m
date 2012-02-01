function s = squish(a, dims)
% squish Remove singleton dimensions more consistently than squeeze()
%   b = squish(a) returns an array b with the same elements as a but with
%   all the singleton dimensions removed.  Unlike the built-in squeeze()
%   function, squish() will reshape row vectors into column vectors; in
%   this sense it has a more consistent behavior than squeeze().
%
%   b = squish(a, dims) will remove only certain singleton dimensions from
%   a.

if nargin == 1
    dims = size(a) == 1;
end

sizeOfA = size(a);
if any(sizeOfA(dims) > 1)
    error('Some dimensions to squish are not singletons');
end
sizeOfA(dims) = [];
if numel(sizeOfA) == 1
    sizeOfA(2) = 1;
end
s = reshape(a, sizeOfA);