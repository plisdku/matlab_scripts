function B = multTensor(T, A, dim)
% mult(T, A, dim) multiplies tensor T by matrix A over dimension dim
% (A*T)

if dim < 1 || dim > ndims(T)
    error('dim must be valid dimension of T');
end

permutedDims = [dim, 1:dim-1, dim+1:ndims(T)];
sz = size(T);
szPermuted = sz(permutedDims);


tPerm = permute(T, permutedDims); % Put the desired dimension first
prodTensor = reshape(A*tPerm(:,:), [size(A,1), szPermuted(2:end)]);

B = ipermute(prodTensor, permutedDims);

