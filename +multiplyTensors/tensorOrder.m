function szC = tensorOrder(szA, szB, inA, inB, replaceDims)

if nargin < 5
    replaceDims = [];
end

if nargin < 3
    inA = [];
end

if nargin < 4
    inB = inA;
end

if ~isempty(replaceDims)
    szA(inA) = szB(replaceDims);
else
    szA(inA) = [];
end

outB = 1:numel(szB);
outB([inB replaceDims]) = [];

szC = [szA szB(outB)];