function [fieldPhasors, z, zParallel, zNormal] = phasors(obj, omegas);
% Called by solve.m

eps0 = PhysicalConstants.eps0;
mu0 = PhysicalConstants.mu0;
c = PhysicalConstants.c;

z = exp(-1i*omegas*obj.dt);

permittivities = zeros(length(obj.permittivityFuncs), length(omegas));
permeabilities = zeros(length(obj.permeabilityFuncs), length(omegas));
assert(all(size(permittivities) == size(permeabilities)));

for layer = 1:obj.numLayers()
    permittivities(layer,:) = obj.permittivityFuncs{layer}(z);
    permeabilities(layer,:) = obj.permeabilityFuncs{layer}(z);
end

%% z

zNormal = ones(obj.numLayers(),length(omegas));
zFirstLayer = ones(3,length(omegas));

if obj.sourceDirection(3) > 0
    firstLayer = 1;
else
    firstLayer = obj.numLayers();
end

tic
zFirstLayer = calculateWaveZ(omegas, ...
    permittivities(firstLayer,:), ...
    permeabilities(firstLayer,:), ...
    obj.sourceDirection, obj.dxyz, obj.dt);
solveDispersionTime = toc;
fprintf('Time to find first layer k vector: %2.2f s\n', solveDispersionTime);

for nn = 1:obj.numLayers()
    zNormal(nn,:) = calculateNormalZ(z, ...
        permittivities(nn,:), permeabilities(nn,:), ...
        zFirstLayer(1:2,:), obj.dxyz, obj.dt);
end

%% phasors

dz = obj.dxyz(3);
fwdPhase = zNormal(1,:).^((obj.layerOrigins(1)-obj.sourceOrigin(3))/dz);
%bkPhase = zNormal(end,:).^((obj.sourceOrigin(3)-obj.layerOrigins(end))/dz);

xPhase = zFirstLayer(1,:).^(obj.sourceOrigin(1)/obj.dxyz(1));
yPhase = zFirstLayer(2,:).^(obj.sourceOrigin(2)/obj.dxyz(2));
sourcePhasors = obj.makeSourcePhasors();
sourcePhasors = sourcePhasors(:,1:length(omegas));
assert(length(sourcePhasors) == length(omegas));
sourcePhasors = [sourcePhasors(1,:).*fwdPhase.*xPhase.*yPhase; ...
    sourcePhasors(2,:).*fwdPhase.*xPhase.*yPhase; ...
    sourcePhasors(3,:)./fwdPhase./xPhase./yPhase; ...
    sourcePhasors(4,:)./fwdPhase./xPhase./yPhase];

fieldPhasors = zeros(obj.numLayers()*12, length(omegas));

tic
for ff = 1:length(omegas)
    fieldPhasors(:,ff) = obj.solveStack(permittivities(:,ff), ...
        permeabilities(:,ff), sourcePhasors(:,ff), z(ff), ...
        zFirstLayer(1:2,ff), zNormal(:,ff));
end
solveStackTime = toc;
fprintf('Time to calculate multilayer phasors: %2.2f s\n', solveStackTime);

zParallel = zFirstLayer(1:2,:);





