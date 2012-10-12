function B = multTensor(T, A, dim)
% mult(T, A, dim) multiplies tensor T by matrix A over dimension dim
% (A*T)

if dim < 1
    error('dim must be valid dimension of T');
end

if dim > ndims(T)
    if ~isvector(A)
        error('A must be a vector in order to multiply singular dimension');
    end
    
    A = reshape(A, [ones(1, dim-1), numel(A)]);
    B = bsxfun(@times, T, A);
else
    % This is the usual case for the multiplication.
    
    permutedDims = [dim, 1:dim-1, dim+1:ndims(T)];
    sz = size(T);
    szPermuted = sz(permutedDims);


    tPerm = permute(T, permutedDims); % Put the desired dimension first
    prodTensor = reshape(A*tPerm(:,:), [size(A,1), szPermuted(2:end)]);

    B = ipermute(prodTensor, permutedDims);
end