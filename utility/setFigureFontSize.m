function setFigureFontSize(size)

set(gca, 'FontSize', size);
set(get(gca, 'YLabel'), 'FontSize', size);
set(get(gca, 'XLabel'), 'FontSize', size);