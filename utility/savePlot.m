function savePlot(filename);
% savePlot          Save the current figure in jpg, eps and fig formats.
%   Usage: savePlot(filename) 

saveas(gcf, [filename, '.png'], 'png');
saveas(gcf, [filename, '.eps'], 'epsc');
saveas(gcf, [filename, '.fig'], 'fig');

