function d = centeredDiff(A, varargin)
% dA = centeredDiff(A)    Like diff(), but with a centered difference.  The
% difference is taken along the first dimension of A with more than one
% entry.
%
% dA = centeredDiff(A, dim) calculates the centered difference along the
% given dimension.

if nargin == 1
    dim = 1;
    if ndims(A) == 2
        if size(A,1) == 1
            dim = 2;
        end
    end
else
    dim = varargin{1};
end

if size(A,dim) == 1
    d = 0*A;
else
    allDims = 1:ndims(A);
    permutation = [dim, setdiff(allDims, dim)];
    dePermutation = [2:dim, 1, dim+1:ndims(A)];

    A = permute(A, permutation);
    d = zeros(size(A));
    d(1,:) = A(2,:) - A(1,:);
    d(end,:) = A(end,:) - A(end-1,:);
    d(2:end-1,:) = 0.5*(A(3:end,:) - A(1:end-2,:));
    d = permute(d, dePermutation);
end