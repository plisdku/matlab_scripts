function [C szC] = txca(A, ndimA, coordA, B, inA, inB)
% txca   Multiply tensor by cell array of matrices
%
% Usage:
%   C = txca(A, ndimA, coordA, B, inA, inB)
%

import multiplyTensors.*

szA = tsize(A, ndimA);

validateArguments(A, szA, coordA, B, inA, inB);
outB = 3-inB;

szC = szA;
szC(inA) = size(B{1}, outB);

C = zeros([szC 1 1]);

indices = repmat({':'}, [1 ndimA]);
for nn = 1:numel(B)
    indices{coordA} = nn;
    
    szA_slice = szA;
    szA_slice(coordA) = 1;
    C(indices{:}) = txa(A(indices{:}), ndimA, B{nn}, inA, inB);
end



function validateArguments(A, szA, coordA, B, inA, inB)

if ~iscell(B)
    error('B must be 1-D cell array of matrices');
end

if size(A, coordA) ~= numel(B)
    error('size(A, coordA) must be the same as numel(B)');
end

for nn = 2:numel(B)
    if ~isequal(size(B{nn}), size(B{1}))
        error('All arrays in B must be the same size!');
    end
end

if size(B{1}, inB) ~= size(A, inA)
    error('size(A, inA) must be equal to size(B, inB)');
end