%% TM tester!

boundaryZ = [-100 0 200 400];
epsr = [1 1+0.1i 2 3 4];
mur = [1 2 1 2 1];
Mx = [0 0 1 0];
Jy = [0 1 0 0];

%boundaryZ = [0];
%epsr = [1 1];
%mur = [1 1];
%Jx = [0];
%My = [1];

omega = 2*pi/800;
ky = 2*pi/12000;

sourceHxLeft = 0;
sourceHxRight = 0;

zOut = linspace(-1000, 1000, 500);


%%

outStruct = solveTM('boundaryZ', boundaryZ, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceHxLeft', sourceHxLeft, 'sourceHxRight', sourceHxRight, ...
    'Mx', Mx, 'Jy', Jy, ...
    'Hx', zOut, 'Ey', zOut, 'Ez', zOut, ...
    'Hxz', zOut, 'Eyz', zOut, 'Ezz', zOut, ...
    'Hxy', zOut, 'Eyy', zOut, 'Ezy', zOut);

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
    
    Hxz = gradient(outStruct.Hx(ii), zOut(ii));
    Eyz = gradient(outStruct.Ey(ii), zOut(ii));
    Ezz = gradient(outStruct.Ez(ii), zOut(ii));
    
    assert( relErr(Hxz, outStruct.Hxz(ii)) < 1e-2);
    assert( relErr(Eyz, outStruct.Eyz(ii)) < 1e-2);
    if norm(outStruct.Ezz) ~= 0
        assert( relErr(Ezz, outStruct.Ezz(ii)) < 1e-2);
    end
    
end

fprintf('Surface normal gradients test PASSED\n');

%% Check the surface-parallel gradients

if ky ~= 0
    Hxy = 1i*ky*outStruct.Hx;
    Eyy = 1i*ky*outStruct.Ey;
    Ezy = 1i*ky*outStruct.Ez;

    assert(relErr(Hxy, outStruct.Hxy) < 1e-2);
    assert(relErr(Eyy, outStruct.Eyy) < 1e-2);
    if norm(outStruct.Ezy) ~= 0
        assert(relErr(Ezy, outStruct.Ezy) < 1e-2);
    end
else
    assert(norm(outStruct.Hxy) == 0);
    assert(norm(outStruct.Eyy) == 0);
    assert(norm(outStruct.Ezy) == 0);
end

fprintf('Surface parallel gradients test PASSED\n');

%% Check Maxwell's equations inside the layers

for nn = 1:numLayers
    
    ii = (zOut > intervals(nn) & zOut < intervals(nn+1));
    
    % -iw*Bx = -dy(Ez) + dz(Ey)
    assert( relErr(-1i*omega*mur(nn)*outStruct.Hx(ii), ...
        -outStruct.Ezy(ii) + outStruct.Eyz(ii)) < 1e-5 );
    
    % -iw*Dy = dz(Hx)
    assert( relErr(-1i*omega*epsr(nn)*outStruct.Ey(ii), ...
        outStruct.Hxz(ii)) < 1e-5 );
    
    % -iw*Dz = -dy(Hx)
    if norm(outStruct.Hxy(ii)) ~= 0
        assert( relErr(-1i*omega*epsr(nn)*outStruct.Ez(ii), ...
            -outStruct.Hxy(ii)) < 1e-5 );
    end
end

fprintf('Maxwell equations test PASSED\n');

%% Check the boundary conditions.
% For this, I'll put output positions closer to boundaries.

zBdy = sort([boundaryZ - 1e-6, boundaryZ + 1e-6]);

outStruct = solveTM('boundaryZ', boundaryZ, ...
    'epsr', epsr, 'mur', mur, 'omega', omega, 'ky', ky, ...
    'sourceHxLeft', sourceHxLeft, 'sourceHxRight', sourceHxRight, ...
    'Mx', Mx, 'Jy', Jy, ...
    'Hx', zBdy, 'Ey', zBdy, 'Ez', zBdy, ...
    'Hxz', zBdy, 'Eyz', zBdy, 'Ezz', zBdy, ...
    'Hxy', zBdy, 'Eyy', zBdy, 'Ezy', zBdy);

for bb = 1:numBoundaries
    
    iLeft = find(zBdy < boundaryZ(bb), 1, 'last');
    iRight = find(zBdy > boundaryZ(bb), 1, 'first');
    
    deltaHx = diff(outStruct.Hx([iLeft, iRight]));
    deltaEy = diff(outStruct.Ey([iLeft, iRight]));
    
    deltaHE = [deltaHx; deltaEy];
    JM = [Jy(bb); Mx(bb)];
    
    assert( norm(JM + deltaHE) < 0.01 ); % H+ - H- = -J
    
    deltaDz = epsr(bb+1)*outStruct.Ez(iRight) - ...
        epsr(bb)*outStruct.Ez(iLeft);
    
    if Jy(bb) == 0
        assert(norm(deltaDz) < 1e-6); % otherwise it should jump with Ex
    end
    
end

fprintf('Boundary conditions test PASSED\n');



