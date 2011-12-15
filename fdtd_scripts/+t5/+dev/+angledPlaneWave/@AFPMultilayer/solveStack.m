% Single-frequency solver
function phasors = solveStack(o, permittivities, permeabilities, ...
    sourcePhasors, z, zParallel, zNormal)
dx = o.dxyz(1);
dy = o.dxyz(2);
dz = o.dxyz(3);

%% System matrix

A = zeros(8*o.numLayers() + 4*length(o.boundariesE), ...
    12*o.numLayers());
row = 1;

%% Maxwell's equations
% There are only four degrees of freedom in the six field components, per
% layer.  I omit Ampere's law for the tangential E field components.

for nn = 1:o.numLayers()
    eps = permittivities(nn);
    mu = permeabilities(nn);
    
    zFw = [zParallel(1), zParallel(2), zNormal(nn)];
    zBk = [zParallel(1), zParallel(2), 1/zNormal(nn)];
    
    %A(row, [o.eForward(1,nn), o.hForward(3,nn), o.hForward(2,nn)]) = ...
    %    [-(z-1)*eps/dt, (1-1/zFw(2))/dy, -(1-1/zFw(3))/dz];
    %A(row+1, [o.eBackward(1,nn), o.hBackward(3,nn), o.hBackward(2,nn)]) = ...
    %    [-(z-1)*eps/dt, (1-1/zBk(2))/dy, -(1-1/zBk(3))/dz];
    
    %A(row, [o.eForward(2,nn), o.hForward(1,nn), o.hForward(3,nn)]) = ...
    %    [-(1-1/z)*eps/dt, (1-1/zFw(3))/dz, -(1-1/zFw(1))/dx];
    %A(row+1, [o.eBackward(2,nn), o.hBackward(1,nn), o.hBackward(3,nn)]) = ...
    %    [-(1-1/z)*eps/dt, (1-1/zBk(3))/dz, -(1-1/zBk(1))/dx];
    
    A(row, [o.hForward(1,nn), o.eForward(3,nn), o.eForward(2,nn)]) = ...
        [-(1-1/z)*mu/o.dt, -(zFw(2)-1)/dy, (zFw(3)-1)/dz];
    A(row+1, [o.hBackward(1,nn), o.eBackward(3,nn), o.eBackward(2,nn)]) = ...
        [-(1-1/z)*mu/o.dt, -(zBk(2)-1)/dy, (zBk(3)-1)/dz];
    
    A(row+2, [o.hForward(2,nn), o.eForward(1,nn), o.eForward(3,nn)]) = ...
        [-(1-1/z)*mu/o.dt, -(zFw(3)-1)/dz, (zFw(1)-1)/dx];
    A(row+3, [o.hBackward(2,nn), o.eBackward(1,nn), o.eBackward(3,nn)]) = ...
        [-(1-1/z)*mu/o.dt, -(zBk(3)-1)/dz, (zBk(1)-1)/dx];
    
    A(row+4, [o.eForward(3,nn), o.hForward(2,nn), o.hForward(1,nn)]) = ...
        [-(z-1)*eps/o.dt, (1-1/zFw(1))/dx, -(1-1/zFw(2))/dy];
    A(row+5, [o.eBackward(3,nn), o.hBackward(2,nn), o.hBackward(1,nn)]) = ...
        [-(z-1)*eps/o.dt, (1-1/zBk(1))/dx, -(1-1/zBk(2))/dy];
    
    A(row+6, [o.hForward(3,nn), o.eForward(2,nn), o.eForward(1,nn)]) = ...
        [-(1-1/z)*mu/o.dt, -(zFw(1)-1)/dx, (zFw(2)-1)/dy];
    A(row+7, [o.hBackward(3,nn), o.eBackward(2,nn), o.eBackward(1,nn)]) = ...
        [-(1-1/z)*mu/o.dt, -(zBk(1)-1)/dx, (zBk(2)-1)/dy];
    
    row = row + 8;
end

%% Boundary conditions

for nn = 1:length(o.boundariesE)
    thisBoundary = o.boundariesE(nn)-o.layerOrigins(nn);
    nextBoundary = o.boundariesE(nn)-o.layerOrigins(nn+1);
    
    A(row, [o.eForward(1,nn), o.eBackward(1,nn), ...
        o.eForward(1,nn+1), o.eBackward(1,nn+1)]) = ...
        ...
        [zNormal(nn)^(thisBoundary/dz), ...
        zNormal(nn)^(-thisBoundary/dz), ...
        -zNormal(nn+1)^(nextBoundary/dz), ...
        -zNormal(nn+1)^(-nextBoundary/dz)];
    
    A(row+1, [o.eForward(2,nn), o.eBackward(2,nn), ...
        o.eForward(2,nn+1), o.eBackward(2,nn+1)]) = ...
        ...
        [zNormal(nn)^(thisBoundary/dz), ...
        zNormal(nn)^(-thisBoundary/dz), ...
        -zNormal(nn+1)^(nextBoundary/dz), ...
        -zNormal(nn+1)^(-nextBoundary/dz)];
    
    row = row + 2;
end

for nn = 1:length(o.boundariesH)
    thisBoundary = o.boundariesH(nn)-o.layerOrigins(nn);
    nextBoundary = o.boundariesH(nn)-o.layerOrigins(nn+1);
    
    A(row, [o.hForward(1,nn), o.hBackward(1,nn), ...
        o.hForward(1,nn+1), o.hBackward(1,nn+1)]) = ...
        ...
        [zNormal(nn)^(thisBoundary/dz), ...
        zNormal(nn)^(-thisBoundary/dz), ...
        -zNormal(nn+1)^(nextBoundary/dz), ...
        -zNormal(nn+1)^(-nextBoundary/dz)];
    
    A(row+1, [o.hForward(2,nn), o.hBackward(2,nn), ...
        o.hForward(2,nn+1), o.hBackward(2,nn+1)]) = ...
        ...
        [zNormal(nn)^(thisBoundary/dz), ...
        zNormal(nn)^(-thisBoundary/dz), ...
        -zNormal(nn+1)^(nextBoundary/dz), ...
        -zNormal(nn+1)^(-nextBoundary/dz)];
    
    row = row + 2;
end
%% Sources!  Extract some fields to source.

sourcePhasors = reshape(sourcePhasors, [], 1); % make column vector

B = A(:,o.nonSourceIndices);
sourceColumns = A(:,o.sourceIndices);

%figure(3)
%imagesc(log(abs(B)), [0, 30])
%colorbar
%pause

warningState = warning();
warning off;
nonSourcePhasors = B \ -sourceColumns*sourcePhasors;
warning(warningState);

%if any(isnan(nonSourcePhasors))
%    warning('NaN');
%end

phasors = zeros(12*length(permittivities),1);
phasors(o.sourceIndices) = sourcePhasors;
phasors(o.nonSourceIndices) = nonSourcePhasors;

%keyboard
