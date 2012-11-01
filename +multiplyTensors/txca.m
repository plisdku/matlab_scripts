function [C szC] = txca(A, B, inA, coordA)
% C = txca(A, B, inA, coordA)   Multiply tensor by cell array of matrices
%

validateArguments(A, B, inA, coordA);

szC = size(A);
szC(inA) = size(B{1}, 1);

C = zeros([szC 1 1]);

indices = repmat({':'}, [1 ndims(C)]);
for nn = 1:numel(B)
    indices{coordA} = nn;
    C(indices{:}) = multiplyTensors.txa(A(indices{:}), B{nn}, inA);
end



function validateArguments(A, B, inA, coordA)

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

if size(B{1}, 2) ~= size(A, inA)
    error('size(A, inA) must be equal to number of columns of elements of B');
end