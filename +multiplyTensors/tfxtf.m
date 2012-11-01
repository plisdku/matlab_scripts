function [C szC] = tfxtf(A, B, coordA, coordB, inA, inB, replaceDims)
% C = tfxtf(A, B, inA, inB, coordA, coordB)   Tensor product of arrays of
% tensors A and B
%
% coordA and coordB are the indices of "spatial" dimensions of the tensor
% fields.  Tensor products will be taken pointwise over these indices.

if nargin < 5
    inA = [];
    inB = [];
elseif nargin < 6
    inB = inA;
end

if nargin < 7
    replaceDims = [];
end

szA = size(A);
szB = size(B);

validateArguments(szA, szB, coordA, coordB, inA, inB, replaceDims);

% Divide the indices of A and B into coordinate and other indices.
% We'll do a permutation and some looping.

% Forward and inverse permutation.  We seek to put the coordinate indices
% of A and B in front of the tensorial indices.
[pA ipA] = permuteToFront(szA, coordA);
[pB ipB] = permuteToFront(szB, coordB);

% Size of tensor residing at each coordinate of A or B
szAA = szA(setdiff(1:numel(szA), coordA));
szBB = szB(setdiff(1:numel(szB), coordB));

% Partition and reshape A and B for multiplication.

A_partitioned = reshape(permute(A, pA), [prod(szA(coordA)) szAA]);
B_partitioned = reshape(permute(B, pB), [prod(szB(coordB)) szBB]);

inAA = ipermute(inA, pA);
inBB = ipermute(inB, pB);

inAA = ipA(inA) - (numel(coordA)-1);
inBB = ipB(inB) - (numel(coordB)-1);

C = elementwiseProduct(A_partitioned, B_partitioned, inAA, inBB);

% Reshape C to have all the coordinates and tensorial indices listed.

szC = size(C);
szC = [szA(coordA), szC(2:end)];
C = reshape(C, [szC 1 1]);

% Permute C to put the coordinates back where they were before, i.e. where
% they were in A, if possible.
%
% Presently, indices of C are [coordA outA outB].  Put them in original
% order again:

outA = setdiff(1:ndims(A), [coordA inA]);
outB = setdiff(1:ndims(B), [coordB inB]);
[unused, pC] = sort([coordA, outA, outB+ndims(A)]);

C = permute(C, pC);
% Now C is indexed in the same order as A originally, then B, less the
% inner product dimensions.  At last we can handle the replacement
% dimensions.

if ~isempty(replaceDims)
    
    % Calculate newRemainingDims.
    % Among the indices of C right now, these are the ones that may have to
    % interleave with outer product dimensions from B.  The important thing
    % is to get them in the right order.
    orderA = 1:ndims(A);
    orderAtrunc = [setdiff(orderA, inA), inA];
    invA = invPerm(orderAtrunc);
    newRemainingDims = invA(setdiff(orderA, inA));
    
    % Calculate newReplaceDims.
    % B has lost its inner product dimensions and coordinate dimensions.
    % The outer product dimensions must be re-ordered to drop in between
    % dimensions from A.  newReplaceDims is an ordering of the outer
    % product indices of B.
    orderB = [replaceDims, setdiff(1:ndims(B), replaceDims)];
    invB = invPerm(orderB);
    newReplaceDims = invB(replaceDims);
    
    % When we provide replaceDims, each other product index of B replaces
    % an inner product index from A.  Thus the size of C is the same as the
    % size of A.  Create a permutation order that puts the replacement
    % indices of B (its outer product indices in the replaceDims order)
    % into the inner product indices, and places the remaining coordinate
    % and outer product indices of A in the other slots.
    orderA(inA) = ndims(A) - numel(inA) + newReplaceDims;
    orderA(setdiff(orderA, inA)) = newRemainingDims;
    
    C = permute(C, orderA);
end







function validateArguments(szA, szB, coordA, coordB, inA, inB, replaceDims)

if numel(unique([coordA inA])) < numel([coordA inA])
    error('Coordinate and inner product indices of A must be distinct');
end

if numel(unique([coordB inB])) < numel([coordB inB])
    error('Coordinate and inner product indices of B must be distinct');
end

if numel(unique([coordB replaceDims])) < numel([coordB replaceDims])
    error('Coordinate and replacement indices of B must be distinct');
end

if numel(coordA) ~= numel(coordB)
    error('A and B must have the same number of coordinate indices');
end

if any(szA(coordA) ~= szB(coordB))
    error('A and B must have the same size over all coordinate indices');
end


function [p, q] = permuteToFront(sz, coord)

p = [coord setdiff(1:numel(sz), coord)];
q = invPerm(p);


function jj = invPerm(ii)

[qq ord] = sort(ii);

ascending = 1:numel(ii);
jj = ascending(ord);


function C = elementwiseProduct(A, B, inA, inB)

assert(size(A,1) == size(B,1));

outA = setdiff(2:ndims(A), inA);
outB = setdiff(2:ndims(B), inB);

szA = size(A);
szB = size(B);
szC = [szA([1, outA]) szB(outB)];

C = zeros([szC 1 1]); % pad with trailing zeros just in case
colons = repmat({':'}, ndims(C)-1);

for nn = 1:size(A,1)
    
    a = reshape(A(nn,:), [szA(2:end) 1 1]);
    b = reshape(B(nn,:), [szB(2:end) 1 1]);
    
    %C(nn, colons{:}) = txt(a, b, inA-1, inB-1);
    
    C(nn, colons{:}) = multiplyTensors.txt_lowLevel(a, b, inA-1, inB-1, [], ...
        outA-1, outB-1, szA(2:end), szB(2:end), szC(2:end));
    
end