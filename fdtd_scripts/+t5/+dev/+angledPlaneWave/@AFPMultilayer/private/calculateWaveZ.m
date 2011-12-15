function [zVector,zAxial] = calculateWaveZ(omegas, permittivity, ...
	permeability, direction, dxyz, dt)

%permittivity = epsVacuum;
%permeability = mu;
%z = zs;

direction = reshape(direction, [], 1);
dxyz = reshape(dxyz, [], 1);
assert(size(permittivity,2) == length(omegas));
assert(size(permeability,2) == length(omegas));

%z = reshape(z,[],1);
permittivity = reshape(permittivity,[],1);
permeability = reshape(permeability,[],1);
omegas = reshape(omegas,[],1);


% The dispersion relation is
%
% mu*eps*(z-1)*(1-1/z)/dt^2 =
%   (zx-1)(1-1/zx)/dx^2 +
%   (zy-1)(1-1/zy)/dy^2 +
%   (zz-1)(1-1/zz)/dz^2
%

timePart = permittivity.*permeability.*sin(omegas*dt/2).^2./dt^2;
%kPart = @(k0) sum(sin(k0*direction.*dxyz/2).^2./dxyz.^2);

kPart = @(k0) sin(0.5*k0*direction(1)*dxyz(1)).^2/(dxyz(1)^2) + ...
    sin(0.5*k0*direction(2)*dxyz(2)).^2/(dxyz(2)^2) + ...
    sin(0.5*k0*direction(3)*dxyz(3)).^2/(dxyz(3)^2);

%%

dxMean = dot(direction.^2, dxyz);

kLinear = (2/dxMean)*asin(sqrt((dxMean/dt)^2*permittivity.*permeability...
    .*sin(omegas*dt/2).^2));

%%

% determine whether a z0 is lossy or not
lossThreshold = 1e-10;
isLossyGainy = @(z) abs(abs(z) - 1.0) >= lossThreshold;

wavevector = zeros(3,length(omegas));
kGuesses = sqrt(permittivity.*permeability.*omegas.^2);
k = 0*kGuesses;

%%

options = optimset('Display', 'off');
for nn = 1:length(omegas)
    kGuess = kLinear(nn);
    fun = @(k0) timePart(nn) - kPart(k0);
    k0 = fsolve(fun,kGuess,options);
    k(nn) = k0;
    wavevector(:,nn) = direction*k0;
end

%% Calculate z
% Fundamentally there should be a pair of k vectors, thus a pair of triples
% of z (zx, zy and zz), that satisfy the dispersion relation at a given
% angle.  Using arcsin can give me one of four k vectors, but if I
% exponentiate them I should still get a valid z, I think.

zVector = exp(1i*wavevector.*repmat(dxyz, 1, length(omegas)));

% I will require that z represent either decrease in amplitude or increase
% in phase along the forward direction.

gainThreshold = 1e-10;
isGainy = @(z,direction) ...
    abs(z.^repmat(direction,1,length(z))) > 1.0 + gainThreshold;
isBackwardPhase = @(z, direction) ...
    angle(z.^repmat(direction,1,length(z))) < 0;

zVector(isBackwardPhase(zVector,direction)) = ...
    1./zVector(isBackwardPhase(zVector,direction));

zVector(isGainy(zVector,direction)) = ...
    1./zVector(isGainy(zVector,direction));

%{
gainyOnes = isGainy(zVector,direction);
backwardsPhaseOnes = ~isGainy(zVector,direction) & ...
    isBackwardPhase(zVector,direction);

zVector(gainyOnes) = 1./zVector(gainyOnes);
zVector(backwardsPhaseOnes) = 1./zVector(backwardsPhaseOnes);
%}

%% Calculate the corresponding z for kLinear

zAxial = exp(1i*kLinear*dxMean);
zAxial(angle(zAxial) < 0) = 1./zAxial(angle(zAxial)<0);
zAxial(abs(zAxial) > 1.0 + gainThreshold) = ...
    1./zAxial(abs(zAxial) > 1.0 + gainThreshold);

assert(~any(abs(zAxial) > (1.0 + gainThreshold)));

