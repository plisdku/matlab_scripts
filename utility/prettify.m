function prettify(fontName, fontSize)

set(gca, 'FontSize', fontSize);
set(get(gca, 'YLabel'), 'FontSize', fontSize);
set(get(gca, 'XLabel'), 'FontSize', fontSize);


set(gca, 'FontName', fontName);
set(get(gca, 'XLabel'), 'FontName', fontName);
set(get(gca, 'YLabel'), 'FontName', fontName);
set(get(gca, 'Title'), 'FontName', fontName);
