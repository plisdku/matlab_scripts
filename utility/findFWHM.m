function [fwhmLeft, fwhmRight] = findFWHM(data1d, peakCell)
%function [fwhmLeft, fwhmRight] = findFWHM(data1d, peakCell)

leftSide = data1d(1:peakCell);
rightSide = data1d(peakCell:end);

peak = data1d(peakCell);

lcell = find(leftSide < peak/2, 1, 'last');
rcell = find(rightSide < peak/2, 1, 'first');


if (~isempty(rcell))
    fwhmRight = peakCell - 1 + ...
        rcell + (peak/2-rightSide(rcell))/(rightSide(rcell-1)-rightSide(rcell));
    
else
    fwhmRight = peakCell;
end
    
if (~isempty(lcell))
    fwhmLeft = lcell + (peak/2-leftSide(lcell))/(leftSide(lcell+1)-leftSide(lcell));
else
    fwhmLeft = peakCell;
end


