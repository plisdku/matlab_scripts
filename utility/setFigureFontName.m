function setFigureFontName(name)

set(gca, 'FontName', name);
set(get(gca, 'XLabel'), 'FontName', name);
set(get(gca, 'YLabel'), 'FontName', name);
set(get(gca, 'Title'), 'FontName', name);