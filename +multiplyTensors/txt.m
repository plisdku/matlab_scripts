function [C, szC] = txt(A, B, inA, inB, replaceDims)
% C = txt(A, B)          Product of tensors, A ox B
% C = txt(A, B, dims)    A ox B contracted over same dimensions of A and B
% C = txt(A, B, inA, inB) contracted over certain dimensions of A and B
% C = txt(A, B, inA, inB, replaceDims) replacing inA with replaceDims
%

import multiplyTensors.*

if nargin < 3
    inA = [];
    inB = [];
elseif nargin < 4
    inB = inA;
end

if nargin < 5
    replaceDims = [];
end

numDims = max(ndims(A), ndims(B));

szA = ones(1, numDims);
szA(1:ndims(A)) = size(A);
szB = ones(1, numDims);
szB(1:ndims(B)) = size(B);

validateArguments(inA, inB, replaceDims, szA, szB)

% Calculate the tensor order of C.
%
% szC is the tensor-theoretical dimension of C.  If the expected result is
% a scalar, then szC will be an empty array.  If the expected result is a
% vector, then szC will be a number and not actually the size of C
% according to Matlab (scalars in Matlab have ndims == 2).

loseDims = @(sizeVec, dims) sizeVec(setdiff(1:numel(sizeVec), dims));
szC = [loseDims(size(A), inA) loseDims(size(B), inB)];
outA = setdiff(1:numel(szA), inA);
outB = setdiff(1:numel(szB), inB);

[C szC] = txt_lowLevel(A, B, inA, inB, replaceDims, outA, outB, szA, szB, szC);







function validateArguments(inA, inB, replaceDims, szA, szB)

if length(inA) ~= length(inB)
    error('Must contract over same number of dimensions');
end

if ~isequal(szA(inA), szB(inB))
    
    errString = 'Contracting dimensions are not all of equal size!';
    
    for dA = inA
        for dB = inB
            if szA(dA) ~= szB(dB)
                errString = [errString, ...
                    sprintf('\n\tsize(A, %i) = %i and size(B, %i) = %i', ...
                        dA, dB, szA(dA), szB(dB))];
            end
        end
    end
                    
    error(errString);
end

if ~isempty(replaceDims)
    if numel(replaceDims) ~= numel(inB)
        error('numel(inB) must equal numel(replaceDims)');
    end
end

