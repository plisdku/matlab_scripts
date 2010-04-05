function epsilon = drudePermittivityFDTD(lambdas, direction, dxdydzdt,...
    epsinf, omegap, gamma)
% permittivity = drudePermittivityFDTD(lambdas, direction, dxdydzdt,
% epsinf, omegap, gamma) returns the actual permittivity of a simulated
% Drude-model material in FDTD with the given cell discretization given by
% dxdydzdt = [dx dy dz dt] for a wave propagating in the given direction =
% [dirx diry dirz], which can be a unit vector but does not need to be.
% Epsinf is a relative permittivity (units of eps0 = 8.85e-12), omegap is
% the plasmon frequency and gamma is the relaxation constant, gamma =
% 1/tauc.

dx = dxdydzdt(1);
dy = dxdydzdt(2);
dz = dxdydzdt(3);
dt = dxdydzdt(4);

omegas = 2*pi*2.99792458e8./lambdas;
eps0 = 8.85418782e-12;

direction = direction/norm(direction);

prefac = 4/(dt^2);   % this all disappears in linearization
eps = epsinf*eps0;
mu = 4e-7*pi;

epsilon = 0*omegas;

for ii = 1:length(omegas)

    w = omegas(ii);
    dirdx = direction.*[dx dy dz];
    wp = omegap;
    g = gamma;
    j1 = (2-g*dt)/(2+g*dt);
    j2 = (2*eps0*wp^2*dt)/(2+g*dt);
    
    bigConst = j2*dt/2/i/eps*sin(w*dt/2)/(exp(-i*w*dt/2) - j1*exp(i*w*dt/2));
    
    dispLin = @(kr,ki) -w^2 - i*w*wp^2/(g-i*w) + (kr+i*ki).^2/mu/eps;
    dispReal = @(kr,ki) prefac*( ...
        -sin(w*dt/2)^2 + bigConst + ...
        dt^2/mu/eps/dx^2*sin( (kr+i*ki)*dirdx(1)/2).^2 + ...
        dt^2/mu/eps/dy^2*sin( (kr+i*ki)*dirdx(2)/2).^2 + ...
        dt^2/mu/eps/dz^2*sin( (kr+i*ki)*dirdx(3)/2).^2 ...
        );

    ksolve1 = sqrt( (w^2*mu*eps)*(1 - wp^2/(w^2 + i*w*g)) );
    ksolve = [real(ksolve1), imag(ksolve1)];

    kbounds = [0.2*ksolve(1), 10*ksolve(1), 0.4*ksolve(2), 1.2*ksolve(2)];

    minfn = @(x) [real(dispReal(x(1),x(2))),imag(dispReal(x(1),x(2)))];
    %{
    opt = optimset('Display', 'off');
    x = lsqnonlin(minfn, ksolve, kbounds([1 3]), kbounds([2 4]), opt);
    %}
    
    opt = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt');
    x = lsqnonlin(minfn, ksolve,[], [], opt);

    k = x(1) + i*x(2);
    epsilon(ii) = k^2 / (mu*w^2);
end

epsilon = epsilon / eps0;
