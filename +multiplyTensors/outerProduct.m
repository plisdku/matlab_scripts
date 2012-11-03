function [C szC] = outerProduct(A, szA, B, szB)
% outerProduct(A, szA, B, szB)  Calculate tensor outer product A ox B
%
% [C szC] = outerProduct(A, szA, B, szC) will satisfy szC == [szA szB].
%

import multiplyTensors.*

if nargin == 4
    C = reshape(outerProduct(A, szA, B(:)), [szA szB]);
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
