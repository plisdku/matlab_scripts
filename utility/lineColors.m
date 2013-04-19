function [theColors] = lineColors(cmap, numLines)
% lineColors  Assign colors from a color map to lines on plot
%
% Example:
%
% lineColors(jet);
% lineColors(hot);

if nargin == 0
    cmap = varycolor(20);
elseif ischar(cmap)
    cmap = feval(cmap);
end

numColors = size(cmap, 1);

if nargin < 2
    children = findobj(get(gca, 'Children'), 'Type', 'line');
    numLines = numel(children);
end

colorParameter = linspace(0, 1, numLines);

theColors = zeros(numLines, 3);

for ll = 1:numLines
    newColor = spline(linspace(0, 1, numColors), cmap', colorParameter(ll))';
    newColor = max(min(newColor, 1), 0);
    if nargin < 2
        set(children(ll), 'Color', newColor);
    end
    theColors(ll,:) = newColor;
end
