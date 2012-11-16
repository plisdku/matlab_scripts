function [C szC] = tfxtf(A, ndimA, coordA, B, ndimB, coordB, inA, inB, replaceDims)
% tfxtf  Calculate product of tensor fields
%
% [C szC] = tfxtf(A, ndimA, coordA, B, ndimB, coordB, inA, inB, replaceB)
%

szA = multiplyTensors.tsize(A, ndimA);
szB = multiplyTensors.tsize(B, ndimB);

if nargin < 9
    replaceDims = [];
end

if nargin < 7
    inA = [];
end

if nargin < 8
    inB = inA;
end

validateArguments(A, szA, coordA, B, szB, coordB, inA, inB, replaceDims);

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

A_partitioned = reshape(permute(A, pA), [prod(szA(coordA)) szAA 1 1]);
B_partitioned = reshape(permute(B, pB), [prod(szB(coordB)) szBB 1 1]);

inAA = ipA(inA) - numel(coordA);
inBB = ipB(inB) - numel(coordB);

C = elementwiseProduct(A_partitioned, szAA, B_partitioned, szBB, ...
    inAA, inBB);

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
[~, pC] = sort([coordA, outA, outB+ndims(A)]);

C = permute(C, pC);
% Now C is indexed in the same order as A originally, then B, less the
% inner product dimensions.  At last we can handle the replacement
% dimensions.

% Now C is the size it "should have been".
if ~isempty(replaceDims)
    
    subArray = @(A,I) A(I);
    % Indices into B(outB) that select only replacement dims.
    repB = subArray(1:(numel(outB) + numel(coordB)), ...
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
    indices(inA) = numel(outA) + numel(coordA) + repB;
    % Now shuffle over the packed-in indices from outA.
    indices(setdiff(1:numel(szA), inA)) = 1:(numel(outA) + numel(coordA));
    
    C = permute(C, indices(1:ndims(C)));
end








function [p, q] = permuteToFront(sz, coord)

p = [coord setdiff(1:numel(sz), coord)];
q = multiplyTensors.invPerm(p);


function C = elementwiseProduct(A, szA, B, szB, inA, inB)
% elementwiseProduct  Multiply two arrays of tensors elementwise
% 
% C = elementwiseProduct(A, szA, B, szB, inA, inB)
%
% A is of nominal size [N szA] and B is of nominal size [N szB].
% Matlab may report fewer dimensions of it strips trailing singletons.
%
% inA and inB are indices into each element of A and B, neglecting the
% leading dimension.

assert(size(A,1) == size(B,1));
numRows = size(A, 1);

outA = setdiff(1:numel(szA), inA);
outB = setdiff(1:numel(szB), inB);

szC = [szA(outA) szB(outB)];

C = zeros([numRows szC 1 1]); % pad with trailing zeros just in case
colons = repmat({':'}, ndims(C)-1);

for nn = 1:numRows
    
    a = reshape(A(nn,:), [szA 1 1]);
    b = reshape(B(nn,:), [szB 1 1]);
    
    C(nn, colons{:}) = multiplyTensors.txt_fast(a, szA, b, szB, inA, inB, [], ...
        outA, outB);
end





function validateArguments(A, szA, coordA, B, szB, coordB, inA, inB, ...
    replaceDims)

% This will take care of most of the validation issues.
multiplyTensors.validateArguments_txt(A, szA, B, szB, inA, inB, replaceDims);

if ~isequal(szA(coordA), szB(coordB))
    error(['Spatial dimensions of A and B are not equal:\n', ...
        '\tszA(coordA) = %s\n\tszB(coordB) = %s'], ...
        num2str(szA(coordA)), num2str(szB(coordB)));
end


if ~isempty(intersect(coordA, inA))
    error(['Spatial indices of A conflict with inner product indices:\n',...
        '\tcoordA = %s\n\tinA = %s'], num2str(coordA), num2str(inA));
end

if ~isempty(intersect(coordB, inB))
    error(['Spatial indices of B conflict with inner product indices:\n',...
        '\tcoordA = %s\n\tinA = %s'], num2str(coordB), num2str(inB));
end

if ~isempty(intersect(coordB, replaceDims))
    error(['Spatial indices of B conflict with replacement indices:\n', ...
        '\tcoordB = %s\n\treplaceDims = %s'], num2str(coordB), ...
        num2str(replaceDims));
end

    

