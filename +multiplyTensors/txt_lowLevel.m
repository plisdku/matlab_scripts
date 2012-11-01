function [C szC] = txt_lowLevel(A, B, inA, inB, replaceDims, outA, outB, szA, szB, szC)

% Now, how to do the products...

if isempty(inA)   % pure outer product!
    C = outerProduct(A, szA, B, szB);
elseif isempty(szC) % contraction over all indices!
    C = dot(A(:), B(:));
else
    % Divide the indices of A and B into inner and outer product indices.

    
    A_out_in = reshape(permute(A, [outA inA]), ...
        prod(szA(outA)), prod(szA(inA)));
    B_in_out = reshape(permute(B, [inB outB]), ...
        prod(szB(inB)), prod(szB(outB)));
    
    % Now A and B are (multi-) indexed like A(I,J) B(J,K), ready for some
    % big matrix operations over J.  I just need to line up the indices
    % appropriately.
    
    if ndims(szC) >= 2
        C = reshape(A_out_in * B_in_out, szC);
    else
        C = reshape(A_out_in * B_in_out, [szC 1]);
    end
    
    assert(isequal(size(C), szC(1:ndims(C))));
    

    % Lastly, handle any desired permutations.  The idea is that we can
    % allow indices of B to replace inner product indices if we want to.

    if ~isempty(replaceDims)
        
        % Fill some of B's outer product dims into the holes left by A's
        % inner product dims.  Some of B's outer product dims remain.
        % Figure out which ones:
        
        subArray = @(A,I) A(I);
        
        outReplaceDims = subArray(1:numel(outB), ...
            subArray((1:numel(szB)) == replaceDims, outB));
        
        newIndices = 1:numel(szA);
        newIndices(inA) = numel(outA) + outReplaceDims;
        newIndices(outA) = 1:numel(outA);
        newIndices = [newIndices, setdiff(1:numel(szC), newIndices)];
        
        C = permute(C, newIndices(1:ndims(C)));
        szC = szC(newIndices);
    end
    
end


function C = outerProduct(A, szA, B, szB)
% outerProduct(A, szA, b)  Calculate outer product tensor A with vector b
% The "real" size of A, szA, must be specified (to include trailing
% singleton dimensions).
%
% C = outerProduct(A, b) will have ndims(C) == ndims(A) + 1
%
% C = outerProduct(A, szA, B, szB) will calculate the outer product of two
% tensors.
% 

if nargin == 4
    C = reshape(outerProduct(A, szA, B(:)), [szA szB]);
else

    assert(isvector(szA));
    assert(isvector(B));

    C = zeros([szA numel(B)]);

    indicesA = repmat({':'}, 1, numel(szA));
    indicesC = indicesA;
    lastIndex = length(indicesC) + 1;

    for db = 1:length(B)

        indicesC{lastIndex} = db;
        C(indicesC{:}) = A*B(db);
        
    end

end
