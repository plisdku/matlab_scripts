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

for jj = 1:length(x)
    minIndex = max([1, jj-nhood(jj)]);
    maxIndex = min([length(x), jj+nhood(jj)]);
    
    maxVal = max(x(minIndex:maxIndex));
    
    
    if maxVal == x(jj)
        indices = [indices, jj];
    end
    
    %if maxVal == x(jj) && sum(x(minIndex:maxIndex) == maxVal) == 1
    %    indices = [indices, jj];
    %end
    
    %{
    [valLeft, indLeft] = max(x(minIndex:jj));
    [valRight, indRight] = max(x(jj:maxIndex));
    
    %if (val == x(jj))
    if (indLeft == jj-minIndex+1) && (indRight == 1)
    %if jj == minIndex + iii - 1
        indices = [indices, jj];
    end
    %}
end

