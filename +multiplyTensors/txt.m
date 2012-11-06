function [C, szC] = txt(A, ndimA, B, ndimB, inA, inB, replaceDims)
% C = txt(A, szA, B, szB, inA, inB, replaceDims)

import multiplyTensors.*

szA = tsize(A, ndimA);
szB = tsize(B, ndimB);

if nargin < 5
    inA = [];
    inB = [];
elseif nargin < 6
    inB = inA;
end

if nargin < 7
    replaceDims = [];
end

validateArguments_txt(A, szA, B, szB, inA, inB, replaceDims)

szC = tensorOrder(szA, szB, inA, inB, replaceDims);

% Obtain the outer indices.
outA = 1:numel(szA); outA(inA) = [];
outB = 1:numel(szB); outB(inB) = [];

C = txt_fast(A, szA, B, szB, inA, inB, replaceDims, outA, outB);

% This is just a sanity check
if numel(szC) >= 2
    szC_actual = size(C);
    assert(isequal(szC_actual, szC(1:numel(szC_actual))));
end

assert(numel(C) == prod(szC));
