function [zNormal, otherRoots] = calculateNormalZ( ...
    z, permittivity, permeability, zParallel, dxyz, dt)

assert(size(permittivity,2) == length(z));
assert(size(permeability,2) == length(z));
assert(size(zParallel,2) == length(z));

%z = exp(1i*linspace(0,0.5));
z = reshape(z,1,[]);
permittivity = reshape(permittivity,1,[]);
permeability = reshape(permeability,1,[]);
dxyz = reshape(dxyz, [], 1);


% The dispersion relation is
%
% (zPerp-1)(1-1/zPerp)/dx^2 = -(zPar-1)(1-1/zPar)/dx^2 + 
%       mu*eps*(z-1)(1-1/z)/dt^2
%
% Express as a root-finding exercise:
%
% rootsOf( c2*zPerp^2 + c1*zPerp + c0 )
%

c2 = 1;
c1 = -2 - dxyz(3)^2*( ...
    sum(-(zParallel-1).*(1-1./zParallel) ./ ...
        repmat(dxyz(1:2).^2,1,length(z)), 1) + ...
    permittivity.*permeability.*(z-1).*(1-1./z)/dt^2);
c0 = 1;

% Determine whether a z is lossy or not
lossThreshold = 1e-10;

isLossyGainy = @(z) abs(abs(z) - 1.0) >= lossThreshold;

%isLossyGainy2 = @(z1, z2) (abs(z1)-1)*(abs(z2)-1) <= lossThreshold^2;

zNormal = [0*z];
otherRoots = zNormal;
for nn = 1:length(z)
    try
        % Find the roots.  Make sure to report them in a consistent order.
        zRoots = roots([c2 c1(nn) c0]);
        
        % How to do this:
        % 1.  Is the material lossy?  Then it's a no-brainer.  Pick the
        %     root that represents loss in the positive direction.
        % 2.  Otherwise, pick the root that represents positive phase
        %     accumulation in the forward direction.
        
        if isLossyGainy(zRoots(1))
            %'lossygainy'
            if abs(zRoots(1)) < abs(zRoots(2))
                zNormal(:,nn) = zRoots(1);
                otherRoots(:,nn) = zRoots(2);
            else
                zNormal(:,nn) = zRoots(2);
                otherRoots(:,nn) = zRoots(1);
            end
            %zNormal(:,nn) = 0;
        else
            %'notlossygainy'
            if imag(zRoots(1)) > imag(zRoots(2))
                zNormal(:,nn) = zRoots(1);
                otherRoots(:,nn) = zRoots(2);
            else
                zNormal(:,nn) = zRoots(2);
                otherRoots(:,nn) = zRoots(1);
            end
        end
        %zNormal(:,nn) = zRoots(round(rand)+1);
    catch
        warning('Root-finding failed');
    end
end
