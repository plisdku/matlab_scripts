function h = pcolor_centered(varargin)

h = pcolor(varargin{:});
clims = get(gca, 'CLim');
clims = [-max(abs(clims)), max(abs(clims))];
set(gca, 'Clim', clims);