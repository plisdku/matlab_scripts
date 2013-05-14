%% TE tester!

import tmm2.*

boundaryX = [-100 0 200 400];
epsr = [1 1+0.1i 2 3 4];
mur = [1 2 1 2 1];
Jz = [0 0 1 0];
Mx = [0 1 0 0];

%boundaryZ = [0];
%epsr = [1 1];
%mur = [1 1];
%Jx = [0];
%My = [1];

omega = 2*pi/800;
ky = 2*pi/12000;

sourceEzLeft = 0;
sourceEzRight = 0;

xOut = linspace(-1000, 1000, 500);


%%

outStruct = solveTE('boundaryX', boundaryX, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceEzLeft', sourceEzLeft, 'sourceEzRight', sourceEzRight, ...
    'Jz', Jz, 'Mx', Mx, ...
    'Ez', xOut, 'Hx', xOut, 'Hy', xOut, ...
    'Ezy', xOut, 'Hxy', xOut, 'Hyy', xOut, ...
    'Ezx', xOut, 'Hxx', xOut, 'Hyx', xOut);

%% Set up basic testing stuff

relErr = @(a, b) norm(a-b)/norm(b);

assert(xOut(1) < boundaryX(1));
assert(xOut(end) > boundaryX(end));
intervals = [xOut(1), boundaryX, xOut(end)];

numLayers = numel(epsr);
numBoundaries = numel(boundaryX);

%% Check the surface-normal gradients inside the layers

for nn = 1:numLayers
    
    ii = (xOut > intervals(nn) & xOut < intervals(nn+1));
    
    Ezx = gradient(outStruct.Ez(ii), xOut(ii));
    Hxx = gradient(outStruct.Hx(ii), xOut(ii));
    Hyx = gradient(outStruct.Hy(ii), xOut(ii));
    
    assert( relErr(Ezx, outStruct.Ezx(ii)) < 1e-2);
    assert( relErr(Hyx, outStruct.Hyx(ii)) < 1e-2);
    if norm(outStruct.Hxx) ~= 0
        assert( relErr(Hxx, outStruct.Hxx(ii)) < 1e-2);
    end
    
end

fprintf('Surface normal gradients test PASSED\n');

%% Check the surface-parallel gradients

if ky ~= 0
    Ezy = 1i*ky*outStruct.Ez;
    Hxy = 1i*ky*outStruct.Hx;
    Hyy = 1i*ky*outStruct.Hy;

    assert(relErr(Ezy, outStruct.Ezy) < 1e-2);
    assert(relErr(Hxy, outStruct.Hxy) < 1e-2);
    if norm(outStruct.Hyy) ~= 0
        assert(relErr(Hyy, outStruct.Hyy) < 1e-2);
    end
else
    assert(norm(outStruct.Ezy) == 0);
    assert(norm(outStruct.Hxy) == 0);
    assert(norm(outStruct.Hyy) == 0);
end

fprintf('Surface parallel gradients test PASSED\n');

%% Check Maxwell's equations inside the layers

for nn = 1:numLayers
    
    ii = (xOut > intervals(nn) & xOut < intervals(nn+1));
    
    % -iw*Dz = dx(Hy) - dy(Hz)
    assert( relErr(-1i*omega*epsr(nn)*outStruct.Ez(ii), ...
        outStruct.Hyx(ii) - outStruct.Hxy(ii)) < 1e-5 );
    
    % -iw*Bx = -dy(Ez)
    if norm(outStruct.Ezy(ii)) ~= 0
        assert( relErr(-1i*omega*mur(nn)*outStruct.Hx(ii), ...
            -outStruct.Ezy(ii)) < 1e-5 );
    end
    
    % -iw*By = dx(Ez)
    assert( relErr(-1i*omega*mur(nn)*outStruct.Hy(ii), ...
        outStruct.Ezx(ii)) < 1e-5 );
end

fprintf('Maxwell equations test PASSED\n');

%% Check the boundary conditions.
% For this, I'll put output positions closer to boundaries.

xBdy = sort([boundaryX - 1e-6, boundaryX + 1e-6]);

outStruct = tmm2.solveTE('boundaryX', boundaryX, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceEzLeft', sourceEzLeft, 'sourceEzRight', sourceEzRight, ...
    'Jz', Jz, 'Mx', Mx, ...
    'Ez', xBdy, 'Hx', xBdy, 'Hy', xBdy, ...
    'Ezy', xBdy, 'Hxy', xBdy, 'Hyy', xBdy, ...
    'Ezx', xBdy, 'Hxx', xBdy, 'Hyx', xBdy);

for bb = 1:numBoundaries
    
    iLeft = find(xBdy < boundaryX(bb), 1, 'last');
    iRight = find(xBdy > boundaryX(bb), 1, 'first');
    
    deltaEz = diff(outStruct.Ez([iLeft, iRight]));
    deltaHy = diff(outStruct.Hy([iLeft, iRight]));
    
    deltaEH = [deltaEz; deltaHy];
    MJ = [Mx(bb); Jz(bb)];
    
    assert( norm(MJ - deltaEH) < 0.01 );
    
    deltaBx = mur(bb+1)*outStruct.Hx(iRight) - ...
        mur(bb)*outStruct.Hx(iLeft);
    
    if Mx(bb) == 0
        assert(norm(deltaBx) < 1e-6); % otherwise it should jump with Ex
    end
    
end

fprintf('Boundary conditions test PASSED\n');



