%% Helper function

assertClose = @(a,b) assert(norm(a-b) < 1e-9);

%% Easy tests with unit simplex

simplex = [0 0 0; 1 0 0; 0 1 0]'; % three points as column vectors

assertClose(ll.inherit.barycentricCoords(simplex, [0 0 0]'), [1 0 0]');
assertClose(ll.inherit.barycentricCoords(simplex, [1 0 0]'), [0 1 0]');
assertClose(ll.inherit.barycentricCoords(simplex, [0 1 0]'), [0 0 1]');


%% Use a weirder triangle

simplex = [-1 0 0; 2 0 0; 0 1 1]';

assertClose(ll.inherit.barycentricCoords(simplex, [0 0 0]'), [2/3 1/3 0]');
assertClose(ll.inherit.barycentricCoords(simplex, [0 1 1]'), [0 0 1]');

%% Multiple points at once

assertClose(ll.inherit.barycentricCoords(simplex, [0 0 0; 0 1 1]'), [2/3 1/3 0; 0 0 1]');


