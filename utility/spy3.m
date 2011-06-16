function spy3(matrix, varargin)
% spy3(matrix, varargin)

[ii, jj, vv] = find(matrix);

if nargin > 1
    plot3(jj,ii,vv, varargin{:});
else
    plot3(jj,ii,vv,'.', 'MarkerSize', 1);
end

axis ij
view(2)
xlabel(sprintf('nz = %i', length(vv)));
xlim([0, size(matrix,2)]);
ylim([0, size(matrix, 1)]);