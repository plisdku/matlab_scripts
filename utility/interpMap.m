function cmap = interpMap(oldColorMap, numSteps)
% cmap = interpMap(oldColorMap, numSteps)
%
% Usage: colormap(interpMap(jet, 1000));

cmap = spline(1:size(oldColorMap,1), oldColorMap', ...
    linspace(1,size(oldColorMap,1), numSteps))';

cmap(:) = max(cmap(:), 0);
cmap(:) = min(cmap(:), 1);

