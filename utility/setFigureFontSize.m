function setFigureFontSize(size)

children = get(gcf, 'Children');

for cc = 1:length(children)
    if strcmpi(get(children(cc), 'Type'), 'axes')

        ax = children(cc);
        set(ax, 'FontSize', size);
        set(get(ax, 'YLabel'), 'FontSize', size);
        set(get(ax, 'XLabel'), 'FontSize', size);
        
    end
end