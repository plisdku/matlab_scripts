function jj = invPerm(ii)
% invPerm   Invert permutation

[~, ord] = sort(ii);

ascending = 1:numel(ii);
jj = ascending(ord);