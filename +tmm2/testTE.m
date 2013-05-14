%% TE tester!

boundaryZ = [-100 0 200 400];
epsr = [1 1+0.1i 2 3 4];
mur = [1 2 1 2 1];
Jx = [0 0 1 0];
My = [0 1 0 0];

%boundaryZ = [0];
%epsr = [1 1];
%mur = [1 1];
%Jx = [0];
%My = [1];

omega = 2*pi/800;
ky = 2*pi/12000;

sourceExLeft = 0;
sourceExRight = 0;

zOut = linspace(-1000, 1000, 500);


%%

outStruct = solveTE('boundaryZ', boundaryZ, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceExLeft', sourceExLeft, 'sourceExRight', sourceExRight, ...
    'Jx', Jx, 'My', My, ...
    'Ex', zOut, 'Hy', zOut, 'Hz', zOut, ...
    'Exz', zOut, 'Hyz', zOut, 'Hzz', zOut, ...
    'Exy', zOut, 'Hyy', zOut, 'Hzy', zOut);

%% Set up basic testing stuff

relErr = @(a, b) norm(a-b)/norm(b);

assert(zOut(1) < boundaryZ(1));
assert(zOut(end) > boundaryZ(end));
intervals = [zOut(1), boundaryZ, zOut(end)];

numLayers = numel(epsr);
numBoundaries = numel(boundaryZ);

%% Check the surface-normal gradients inside the layers

for nn = 1:numLayers
    
    ii = (zOut > intervals(nn) & zOut < intervals(nn+1));
    
    Exz = gradient(outStruct.Ex(ii), zOut(ii));
    Hyz = gradient(outStruct.Hy(ii), zOut(ii));
    Hzz = gradient(outStruct.Hz(ii), zOut(ii));
    
    assert( relErr(Exz, outStruct.Exz(ii)) < 1e-2);
    assert( relErr(Hyz, outStruct.Hyz(ii)) < 1e-2);
    if norm(outStruct.Hzz) ~= 0
        assert( relErr(Hzz, outStruct.Hzz(ii)) < 1e-2);
    end
    
end

fprintf('Surface normal gradients test PASSED\n');

%% Check the surface-parallel gradients

if ky ~= 0
    Exy = 1i*ky*outStruct.Ex;
    Hyy = 1i*ky*outStruct.Hy;
    Hzy = 1i*ky*outStruct.Hz;

    assert(relErr(Exy, outStruct.Exy) < 1e-2);
    assert(relErr(Hyy, outStruct.Hyy) < 1e-2);
    if norm(outStruct.Hzy) ~= 0
        assert(relErr(Hzy, outStruct.Hzy) < 1e-2);
    end
else
    assert(norm(outStruct.Exy) == 0);
    assert(norm(outStruct.Hyy) == 0);
    assert(norm(outStruct.Hzy) == 0);
end

fprintf('Surface parallel gradients test PASSED\n');

%% Check Maxwell's equations inside the layers

for nn = 1:numLayers
    
    ii = (zOut > intervals(nn) & zOut < intervals(nn+1));
    
    % -iw*Dx = dy(Hz) - dz(Hy)
    assert( relErr(-1i*omega*epsr(nn)*outStruct.Ex(ii), ...
        outStruct.Hzy(ii) - outStruct.Hyz(ii)) < 1e-5 );
    
    % -iw*By = -dz(Ex)
    assert( relErr(-1i*omega*mur(nn)*outStruct.Hy(ii), ...
        -outStruct.Exz(ii)) < 1e-5 );
    
    % -iw*Bz = dy(Ex)
    if norm(outStruct.Exy(ii)) ~= 0
        assert( relErr(-1i*omega*mur(nn)*outStruct.Hz(ii), ...
            outStruct.Exy(ii)) < 1e-5 );
    end
end

fprintf('Maxwell equations test PASSED\n');

%% Check the boundary conditions.
% For this, I'll put output positions closer to boundaries.

zBdy = sort([boundaryZ - 1e-6, boundaryZ + 1e-6]);

outStruct = solveTE('boundaryZ', boundaryZ, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceExLeft', sourceExLeft, 'sourceExRight', sourceExRight, ...
    'Jx', Jx, 'My', My, ...
    'Ex', zBdy, 'Hy', zBdy, 'Hz', zBdy, ...
    'Exz', zBdy, 'Hyz', zBdy, 'Hzz', zBdy, ...
    'Exy', zBdy, 'Hyy', zBdy, 'Hzy', zBdy);

for bb = 1:numBoundaries
    
    iLeft = find(zBdy < boundaryZ(bb), 1, 'last');
    iRight = find(zBdy > boundaryZ(bb), 1, 'first');
    
    deltaEx = diff(outStruct.Ex([iLeft, iRight]));
    deltaHy = diff(outStruct.Hy([iLeft, iRight]));
    
    deltaEH = [deltaEx; deltaHy];
    MJ = [My(bb); Jx(bb)];
    
    assert( norm(MJ - deltaEH) < 0.01 );
    
    deltaBz = mur(bb+1)*outStruct.Hz(iRight) - ...
        mur(bb)*outStruct.Hz(iLeft);
    
    if My(bb) == 0
        assert(norm(deltaBz) < 1e-6); % otherwise it should jump with Ex
    end
    
end

fprintf('Boundary conditions test PASSED\n');



