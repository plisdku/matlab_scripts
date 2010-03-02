% test z drude

lambdas = 300e-9:10e-9:1200e-9;
omegas = 2*pi*3e8./lambdas;

eps0 = 8.85e-12;
mu0 = 4e-7*pi;

dx = 5e-9;
dt = 0.99*dx/sqrt(3)/3e8;

omegap = 4.0217e15;
tauc = 9.0441e-15;
epsinf = 12.9898;

[epsPalik, lamPalik] = getPermittivity('Au', [lambdas(1), lambdas(end)]);
epsTheory = drudePermittivity(lambdas, epsinf, omegap, 1/tauc);
[numerator, denominator] = zDrude4(dt, 'omegap', omegap, 'tauc', tauc, ...
    'epsinf', epsinf);
zVals = exp(-i*omegas*dt);
zPermittivity = polyval(numerator, zVals)./polyval(denominator, zVals);

figure(1)
clf
plot(lambdas*1e9, real(epsTheory), lambdas*1e9, imag(epsTheory));
hold on
plot(lamPalik*1e9, real(epsPalik), 'o', lamPalik*1e9, imag(epsPalik), 'o');
plot(lambdas*1e9, real(zPermittivity), 'r.',...
    lambdas*1e9, imag(zPermittivity), 'r.');


figure(2)
clf
plot(lambdas*1e9, real(zPermittivity), lambdas*1e9, imag(zPermittivity));

