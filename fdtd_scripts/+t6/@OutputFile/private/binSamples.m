function [N, bins, indices] = binSamples(y, x)
% [bins, firstIndices] = binSamples(y,x)
%
% put samples y into bins separated by x

% The bins from histc are:
% [x(1), x(2))
% [x(2), x(3))
% ...
% [x(end-1), x(end))
% [x(end), x(end)]
% so the last bin contains samples that are exactly on the right endpoint.

% In my case I am goig to toss the last bin into the second to last bin
% because I really want to bracket each output sample in an interval.

[N, bins] = histc(y,x);
N(end-1) = N(end-1) + N(end);
N = N(1:end-1);
bins(bins == length(x)) = length(x)-1;

indices = cumsum([1 N]);
