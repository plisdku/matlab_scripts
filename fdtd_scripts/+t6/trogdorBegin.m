function trogdorBegin(varargin)
%trogdorBegin Begin Trogdor simulation description
%   trogdorBegin('Bounds', [0 0 0 100 100 100], 'NumCells', [20 20 20], 
%       'Duration', 100, 'Courant', 0.99) declares a simulation where each
%       Yee cell has size [5 5 5] and the physical duration of the
%       simulation is at least 100 (in simulation units).  The Courant
%       parameter adjusts the timestep based on the spatial step.

X.Bounds = [0 0 0 0 0 0];
X.NumCells = [1 1 1];
X.Duration = 1;
X.Courant = 0.99;

X = parseargs(X, varargin{:});

if any(X.NumCells < 1)
    error('All elements of NumCells must be at least 1.');
end

dxyz = (X.Bounds(4:6) - X.Bounds(1:3)) ./ X.NumCells;
dxyz(dxyz == 0) = 1;

assert(numel(dxyz) == 3);

dimensions = [];
for xyz = 1:3
    if X.Bounds(xyz) < X.Bounds(xyz+3) && X.NumCells(xyz) > 1
        dimensions = [dimensions, xyz];
    end
end

dt = X.Courant * t6.courant(dxyz(dimensions));
numTimesteps = ceil(X.Duration / dt);

t6.TrogdorSimulation.clear();
sim = t6.TrogdorSimulation.instance();

sim.Dxyz = dxyz;
sim.Dt = dt;
sim.NumT = numTimesteps;
sim.NumCells = X.NumCells;
sim.Grids = [];
sim.CurrentGrid = [];

assert(numel(sim.Dxyz) == 3);

t6.addGrid('Main', [0, 0, 0, X.NumCells-1], X.Bounds(1:3));


