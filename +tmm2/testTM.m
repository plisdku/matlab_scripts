%% TM tester!

boundaryX = [-100 0 200 400];
epsr = [1 1+0.1i 2 3 4];
mur = [1 2 1 2 1];
Mz = [0 0 1 0];
Jy = [0 1 0 0];

%boundaryZ = [0];
%epsr = [1 1];
%mur = [1 1];
%Jx = [0];
%My = [1];

omega = 2*pi/800;
ky = 2*pi/12000;

sourceHzLeft = 0;
sourceHzRight = 0;

xOut = linspace(-1000, 1000, 500);


%%

outStruct = tmm2.solveTM('boundaryX', boundaryX, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceHzLeft', sourceHzLeft, 'sourceHzRight', sourceHzRight, ...
    'Mz', Mz, 'Jy', Jy, ...
    'Hz', xOut, 'Ex', xOut, 'Ey', xOut, ...
    'Hzy', xOut, 'Exy', xOut, 'Eyy', xOut, ...
    'Hzx', xOut, 'Exx', xOut, 'Eyx', xOut);

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
    
    Hzx = gradient(outStruct.Hz(ii), xOut(ii));
    Exx = gradient(outStruct.Ex(ii), xOut(ii));
    Eyx = gradient(outStruct.Ey(ii), xOut(ii));
    
    assert( relErr(Hzx, outStruct.Hzx(ii)) < 1e-2);
    if norm(outStruct.Exx) ~= 0
        assert( relErr(Exx, outStruct.Exx(ii)) < 1e-2);
    end
    assert( relErr(Eyx, outStruct.Eyx(ii)) < 1e-2);
    
end

fprintf('Surface normal gradients test PASSED\n');

%% Check the surface-parallel gradients

if ky ~= 0
    Hzy = 1i*ky*outStruct.Hz;
    Exy = 1i*ky*outStruct.Ex;
    Eyy = 1i*ky*outStruct.Ey;

    assert(relErr(Hzy, outStruct.Hzy) < 1e-2);
    if norm(outStruct.Exy) ~= 0
        assert(relErr(Exy, outStruct.Exy) < 1e-2);
    end
    assert(relErr(Exy, outStruct.Exy) < 1e-2);
else
    assert(norm(outStruct.Hzy) == 0);
    assert(norm(outStruct.Exy) == 0);
    assert(norm(outStruct.Eyy) == 0);
end

fprintf('Surface parallel gradients test PASSED\n');

%% Check Maxwell's equations inside the layers

for nn = 1:numLayers
    
    ii = (xOut > intervals(nn) & xOut < intervals(nn+1));
    
    % -iw*Bx = -dy(Ez) + dz(Ey)
    assert( relErr(-1i*omega*mur(nn)*outStruct.Hz(ii), ...
        -outStruct.Eyx(ii) + outStruct.Exy(ii)) < 1e-5 );
    
    % -iw*Dy = dz(Hx)
    if norm(outStruct.Hzy(ii)) ~= 0
        assert( relErr(-1i*omega*epsr(nn)*outStruct.Ex(ii), ...
            outStruct.Hzy(ii)) < 1e-5 );
    end
    
    % -iw*Dz = -dy(Hx)
    assert( relErr(-1i*omega*epsr(nn)*outStruct.Ey(ii), ...
        -outStruct.Hzx(ii)) < 1e-5 );
end

fprintf('Maxwell equations test PASSED\n');

%% Check the boundary conditions.
% For this, I'll put output positions closer to boundaries.

xBdy = sort([boundaryX - 1e-6, boundaryX + 1e-6]);

outStruct = tmm2.solveTM('boundaryX', boundaryX, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceHzLeft', sourceHzLeft, 'sourceHzRight', sourceHzRight, ...
    'Mz', Mz, 'Jy', Jy, ...
    'Hz', xBdy, 'Ex', xBdy, 'Ey', xBdy, ...
    'Hzy', xBdy, 'Exy', xBdy, 'Eyy', xBdy, ...
    'Hzx', xBdy, 'Exx', xBdy, 'Eyx', xBdy);

for bb = 1:numBoundaries
    
    iLeft = find(xBdy < boundaryX(bb), 1, 'last');
    iRight = find(xBdy > boundaryX(bb), 1, 'first');
    
    deltaHz = diff(outStruct.Hz([iLeft, iRight]));
    deltaEy = diff(outStruct.Ey([iLeft, iRight]));
    
    deltaHE = [deltaHz; deltaEy];
    JM = [Jy(bb); Mz(bb)];
    
    assert( norm(JM + deltaHE) < 0.01 ); % H+ - H- = -J
    
    deltaDx = epsr(bb+1)*outStruct.Ex(iRight) - ...
        epsr(bb)*outStruct.Ex(iLeft);
    
    if Jy(bb) == 0
        assert(norm(deltaDx) < 1e-6); % otherwise it should jump with Ex
    end
    
end

fprintf('Boundary conditions test PASSED\n');



