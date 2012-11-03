function [C szC] = txa(A, ndimA, B, inA, inB)
% txa  Multiply tensor by array
% 
% [C szC] = txa(A, szA, B, inA, inB)

%szA = multiplyTensors.tsize(A, ndimA);

if nargin < 5
    inB = 2;
end

outB = 3-inB;

[C szC] = multiplyTensors.txt(A, ndimA, B, 2, inA, inB, outB);