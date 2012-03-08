function lineColors(cmap)
% lineColors  Assign colors from a color map to lines on plot
%
% Example:
%
% lineColors(jet);
% lineColors(hot);

numColors = size(cmap, 1);

children = findobj(get(gca, 'Children'), 'Type', 'line');

numLines = numel(children);

colorParameter = linspace(0, 1, numLines);

for ll = 1:numLines
    newColor = spline(linspace(0, 1, numColors), cmap', colorParameter(ll))';
    newColor = max(min(newColor, 1), 0);
    set(children(ll), 'Color', newColor);
end