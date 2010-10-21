function d = centeredDiff(A, varargin)

if nargin == 1
    dim = 1;
else
    dim = varargin{1};
end

allDims = 1:ndims(A);
permutation = [dim, setdiff(allDims, dim)];
dePermutation = [2:dim, 1, dim+1:ndims(A)];

A = permute(A, permutation);
d = zeros(size(A));
d(1,:) = A(2,:) - A(1,:);
d(end,:) = A(end,:) - A(end-1,:);
d(2:end-1,:) = 0.5*(A(3:end,:) - A(1:end-2,:));
d = permute(d, dePermutation);
