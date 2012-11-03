function [C szC] = txa(A, szA, B, inA, inB)
% txa  Multiply tensor by array
% 
% [C szC] = txa(A, szA, B, inA, inB)

if nargin < 5
    inB = 2;
end

outB = 3-inB;

[C szC] = multiplyTensors.txt(A, szA, B, size(B), inA, inB, outB);