function indices = findLocalMin(x, neighborhood)
% findLocalMin          find points which are lower than their neighbors.
%   indices = findLocalMin(x, neighborhood) returns indices ii into x of
%   data in vector x which are the min of the subset
%
%       x(ii-neighborhood:ii+neighborhood)
%
%   Neighborhood should be a scalar or a vector the same size as x.

indices = findLocalMax(-x, neighborhood);

%{
if (length(neighborhood) == length(x))
    nhood = neighborhood;
elseif (length(neighborhood) == 1)
    nhood = neighborhood * ones(size(x));
else
    error('Input neighborhood must be scalar or vector of same size as x.');
end

for (jj = 1:length(x))
    minIndex = max([1, jj-nhood(jj)]);
    maxIndex = min([length(x), jj+nhood(jj)]);
    
    
    [val, iii] = min(x(minIndex:maxIndex));
    
    if (val == x(jj))
        indices = [indices, jj];
    end
end
%}