%% Meshgrid orientation
% This demonstrates that meshgrid returns matrices that increase in the
% row direction and increase in the column direction, respectively.

Nrows = 30;
Ncols = 20;

colIndices = repmat(1:Ncols, Nrows, 1);
rowIndices = repmat( (1:Nrows)', 1, Ncols);
[xx, yy] = meshgrid(1:Nrows, 1:Ncols);

% Now xx(1,:) = 1:Nrows and yy(:,1) = (1:Ncols)'.

%% Quiver orientation
figure(1);
quiver(xx, yy, xx, 0*yy);
axis image
xlabel('x = increasing column index');
ylabel('y = increasing row index');
title('quiver(xx, yy, xx, 0*yy)');

figure(2)
quiver(xx, yy, 0*xx, yy);
axis image
xlabel('x = increasing column index');
ylabel('y = increasing row index');
title('quiver(xx, yy, 0*xx, yy)');

%%