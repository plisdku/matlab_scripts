function h = imagesc_centered(varargin)

h = imagesc(varargin{:});
clims = get(gca, 'CLim');
clims = [-max(abs(clims)), max(abs(clims))];
set(gca, 'Clim', clims);