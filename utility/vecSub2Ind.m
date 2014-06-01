function ii = vecSub2Ind(sz, ijk)
% idx = vecSub2Ind(sz, ijk)  Convert subscript array to linear array index
%
% sz is the size of the array to index into, length(sz) == numdims
% ijk is size [N numdims]
% ii is size [N 1].
%
% This is like sub2ind.  If ii, jj, kk are column vectors, then the
% following will give the same result:
%
% sub2ind(sz, ii, jj, kk)
% vecSub2Ind(sz, [ii jj kk]).

numdims = length(sz);

assert(size(ijk,2) == numdims);

sliceSize = [1 cumprod(sz(1:end-1))];

%ii = 1 + (ijk-1)*sliceSize';
ii = 1 + ijk*sliceSize' - sum(sliceSize); % subtle, faster!