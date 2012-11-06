function C = txt_fast(A, szA, B, szB, inA, inB, replaceDims, outA, outB)

import multiplyTensors.*

%if isscalar(A) || isscalar(B)
%    C = A*B;
%    return
%end

% Neither A nor B is a scalar.  
%
% [2 1] representing [2 1] needs to (?) push over to [1 2] (shift -1)
% [2 1] representing [2 1 1] needs to push over to [1 1 2] (shift -1)
% [1 2] representing [1 2 1] needs to push over to [1 1 2] (shift -1)
%
% I need to count the number of trailing 1s in szA and szB to figure out
% the right shift.
%
% How do I undo this at the end?  Let's see...
%
% [3 2 1] ox [4 4 1] should be [3 2 1 4 4 1]
% [1 3 2] ox [1 4 4] should be [1 3 2 1 4 4]
% [1 3 2] ox [1] should be [1 3 2 1] but drops to [1 3 2]
% 

onesA = numTrailingOnes(szA);
onesB = numTrailingOnes(szB);

if onesA > 0
    A = shiftdim(A, -onesA);
    inA_ = wrapIndices(numel(szA), inA + onesA);
    outA_ = wrapIndices(numel(szA), outA + onesA);
    szA_ = szA(wrapIndices(numel(szA), (1:numel(szA)) - onesA));
    assert(szA_(1) == size(A,1));
    orderA = [numel(szA) - onesA + (1:onesA), 1:(numel(szA) - onesA)];
else
    [inA_ outA_ szA_] = deal(inA, outA, szA);
    orderA = 1:numel(szA);
end

if onesB > 0 && ~ismatrix(B) % matrices are always 2D, no probs here
    B = shiftdim(B, -onesB);
    inB_ = wrapIndices(numel(szB), inB + onesB);
    outB_ = wrapIndices(numel(szB), outB + onesB);
    %replaceDims_ = wrapIndices(numel(szB), replaceDims + onesB);
    szB_ = szB(wrapIndices(numel(szB), (1:numel(szB)) - onesB));
    orderB = [numel(szB) - onesB + (1:onesB), 1:(numel(szB) - onesB)];
else
    [inB_ outB_ szB_] = deal(inB, outB, szB);
    orderB = 1:numel(szB);
end
[~, orderA] = sort(orderA(outA_));
[~, orderB] = sort(orderB(outB_));
reorderC = invPerm([orderA, orderB + numel(orderA)]);

%assert(ndims(A) == numel(szA));
%assert(ndims(B) == numel(szB));


if isempty(inA)
    C = permute(outerProduct(A, ndims(A), B, ndims(B)), reorderC);
    
else
    C = txt_faster(A, szA_, B, szB_, inA_, outA_, inB_, outB_);
    
    % Un-permute C from the leading singleton ordering.
    %
    % C is now ordered [ones OA ones OB], though actually all of the ones
    % are inner product indices.  We can cull the orderings now.
    
    
    % Now C is the size it "should have been".
    if ~isempty(replaceDims)
        
        % Indices into B(outB) that select only replacement dims.
        repB = subArray(1:numel(outB), ...
            subArray( (1:numel(szB)) == replaceDims, outB));
        % The effect here is that each element of repB is smaller than
        % corresponding elements of replaceDims by the number of inner
        % product indices of B that were beneath it and got contracted out.
        
        % Now its complement: indices into B(outB) that select only
        % non-replacement dims.
        nonRepB = setdiff(1:numel(outB), repB);
        assert(numel(intersect(repB, nonRepB)) == 0);
        
        % Build the new endex ordering.  The array built here will re-order
        % the dimensions of C_out_out.
        indices = 1:numel(szA);
        indices(inA) = numel(outA) + repB;
        % Now shuffle over the packed-in indices from outA.
        indices(outA) = 1:numel(outA);
        
        %C = permute(C, indices(1:ndims(C)));
        
        % I'm going to try to compose these permutations.
        
        reorderC = reorderC(indices);
    end
    
    if ~isempty(reorderC)
        C = permute(C, reorderC);
    end
    
end


function C = txt_faster(A, szA, B, szB, inA, outA, inB, outB)

% Divide the indices of A and B into inner and outer product indices.
% That is, reorder indices to [outA inA] and [inB outB] for matrix
% multiplication.
A_out_in = reshape(permute(A, [outA inA]), ...
    prod(szA(outA)), prod(szA(inA)));
B_in_out = reshape(permute(B, [inB outB]), ...
    prod(szB(inB)), prod(szB(outB)));

% Here's C, sorted with indices of A before indices of B.  The reshape
% function tacks on extra singletons just in case szC was less than
% two elements long, since Matlab can't reshape an array to less than
% two dimensions.

szC_current = [szA(outA) szB(outB)];
C = reshape(A_out_in*B_in_out, [szC_current 1 1]);



function A = subArray(B, I)
A = B(I);

function n = numTrailingOnes(sz)

n = numel(sz) - find(sz > 1, 1, 'last');


function ii = wrapIndices(numDims, indices)

ii = 1 + mod(indices-1, numDims);
