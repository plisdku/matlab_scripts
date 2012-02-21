%% Test Trogdor grid stuff

import t6.*

trogdorBegin('Bounds', [0 0 0 10 10 0], 'NumCells', [10 10 1], ...
    'PML', [1 1 0 2 2 0]);

sim = simulation();

expectedOuterBounds = [-1 -1 0 12 12 0];
expectedInnerBounds = [0 0 0 10 10 0];
expectedNumCells = [13 13 1];

assert(all(expectedOuterBounds == sim.OuterBounds));
assert(all(expectedInnerBounds == sim.NonPMLBounds));
assert(all(expectedNumCells == sim.NumCells));

%currentGrid = grid();

trogdor_end('XML', 'testParams.xml')
delete('testParams.xml');

fprintf('Test PML bounds: PASSED\n');

