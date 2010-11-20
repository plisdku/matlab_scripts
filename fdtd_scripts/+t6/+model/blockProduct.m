function [C, blockA, blockB] = blockProduct(A, B, blockSizeA, blockSizeB)
%blockProduct    Blockwise product of matrices (poor man's tensor product)
% A block matrix is a matrix comprised of an integral tiling of matrices of
% uniform size, for instance
%
% A = [A1 A2 A3]                       blockwise row vector
% B = [B1; B2; B3]                     blockwise column vector
% C = [C11 C12; C21 C22; C31 C32]      blockwise matrix
% D = [D1]                             blockwise scalar
%
% If the blocks of matrices A, B and C are the same size, blockwise
% multiplication is intuitively defined:
%
% A*D = [A1*D A2*D A3*D];
%
% A*B = [ A1*B1 A2*B1 A3*B1;
%         A1*B2 A2*B2 A3*B2;
%         A1*B3 A2*B3 A3*B3; ]
%
% C*D = [ C11*D  C12*D
%         C21*D  C22*D
%         C31*D  C32*D ]
%
% When two matrices have multiple blocks along a shared dimension, their
% product is not defined here.
%
% To remove ambiguities about matrix size, the block dimensions of each
% matrix must be provided explicitly.
%
% Usage:
%
% C = blockProduct(A, B, blockSizeA, blockSizeB)
%
% Example:
%   A = [eye(2) 2*eye(2)];
%   B = [3*eye(2); 4*eye(2)];
%   blockProduct(A, B, [2 2], [2 2]);
% 
% ans =
% 
%      3     0     6     0
%      0     3     0     6
%      4     0     8     0
%      0     4     0     8
%
% Author: Paul C. Hansen
% pch@stanford.edu

numBlocksA = size(A)./blockSizeA;
numBlocksB = size(B)./blockSizeB;

% Basic error checking: if the user claims to have sent in a block matrix,
% make sure each block is the same integral size; make sure as well that
% the blocks are the right size to be multiplied together.
if blockSizeA(2) ~= blockSizeB(1)
    error('Blocks of matrix A and matrix B do not have the same inner dimensions');
end

if any(round(numBlocksA) ~= numBlocksA)
    error('Matrix A has noninteger number of blocks');
end

if any(round(numBlocksB) ~= numBlocksB)
    error('Matrix B has noninteger number of blocks');
end

if numBlocksA(1) > 1 && numBlocksB(1) > 1
    error('Only one of the input matrices may have multiple row blocks.');
end

if numBlocksA(2) > 1 && numBlocksB(2) > 1
    error('Only one of the input matrices may have multiple column blocks.');
end

% Make some anonymous functions to extract block(m,n) of the input
% matrices.
if numBlocksA(1) > 1
    if numBlocksA(2) > 1
        blockA = @(m,n) A( (1+(m-1)*blockSizeA(1)):m*blockSizeA(1), ...
            (1+(n-1)*blockSizeA(2)):n*blockSizeA(2));
    else
        blockA = @(m,n) A((1+(m-1)*blockSizeA(1)):m*blockSizeA(1), :);
    end
elseif numBlocksA(2) > 1
    blockA = @(m,n) A(:, (1+(n-1)*blockSizeA(2)):n*blockSizeA(2));
else
    blockA = @(m,n) A;
end

if numBlocksB(1) > 1
    if numBlocksB(2) > 1
        blockB = @(m,n) B( (1+(m-1)*blockSizeB(1)):m*blockSizeB(1), ...
            (1+(n-1)*blockSizeB(2)):n*blockSizeB(2));
    else
        blockB = @(m,n) B((1+(m-1)*blockSizeB(1)):m*blockSizeB(1), :);
    end
elseif numBlocksB(2) > 1
    blockB = @(m,n) B(:, (1+(n-1)*blockSizeB(2)):n*blockSizeB(2));
else
    blockB = @(m,n) B;
end

% Calculate the products blockwise.

numBlocks = max(numBlocksA, numBlocksB);
blockSizeC = [blockSizeA(1) blockSizeB(2)];

if issparse(A) || issparse(B)
    sparseSize = numBlocks.*blockSizeC;
    C = sparse(sparseSize(1), sparseSize(2));
else
    C = zeros(numBlocks.*blockSizeC);
end

for row = 1:numBlocks(1)
    for col = 1:numBlocks(2)
        C( (1+(row-1)*blockSizeC(1)):row*blockSizeC(1), ...
            (1+(col-1)*blockSizeC(2)):col*blockSizeC(2) ) = ...
            blockA(row,col) * blockB(row,col);
    end
end



