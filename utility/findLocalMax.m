function indices = findLocalMax(x, neighborhood)
% findLocalMax          find points which are higher than their neighbors.
%   indices = findLocalMax(x, neighborhood) returns indices ii into x of
%   data in vector x which are the max of the subset
%
%       x(ii-neighborhood:ii+neighborhood)
%
%   Neighborhood should be a scalar or a vector the same size as x.

indices = [];

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
    
    
    [val, iii] = max(x(minIndex:maxIndex));
    
    if (val == x(jj))
        indices = [indices, jj];
    end
end

