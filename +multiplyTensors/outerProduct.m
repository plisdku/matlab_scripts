function [C szC] = outerProduct(A, ndimA, B, ndimB)
% outerProduct(A, szA, B, szB)  Calculate tensor outer product A ox B
%
% [C szC] = outerProduct(A, szA, B, szC) will satisfy szC == [szA szB].
%

import multiplyTensors.*

szA = tsize(A, ndimA);

if nargin == 4
    szB = tsize(B, ndimB);
    C = reshape(outerProduct(A, ndimA, B(:)), [szA szB]);
else

    assert(isvector(szA));
    assert(isvector(B));

    C = zeros([szA numel(B)]);

    indicesA = repmat({':'}, 1, numel(szA));
    indicesC = indicesA;
    lastIndex = length(indicesC) + 1;

    for db = 1:length(B)

        indicesC{lastIndex} = db;
        C(indicesC{:}) = A*B(db);
        
    end

end

if nargout > 1
    szC = [szA szB];
end
