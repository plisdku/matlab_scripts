% test zCFSRIPML


dx = 5e-9;
dt = 0.99*dx/sqrt(3)/3e8;

% Parameters for the PML
nPML = 10;  % number of Yee cells of PML
L = nPML*dx;
R0 = 1e-6;  % desired attenuation at normal incidence
R0 = 1e-7;
eps0 = 8.85e-12;
c = 3e8;
m = 3.5;  % polynomial scaling exponent

sigmaMax = -log(R0) * (m+1) * c*eps0/(2*L);
alphaMax = eps0*c/(10*dx);  % use size of scatterer
kappaMax = 5;

% d goes from 0 to 1.
kappaFn = @(d) 1 + (kappaMax-1)*max(d,0).^3;
sigmaFn = @(d) sigmaMax*max(d,0).^3;
alphaFn = @(d) alphaMax + 0*max(d,0);

sFn = @(kappa, sigma, alpha, omega) kappa + sigma ./ (alpha - 1i*omega*eps0);

d = linspace(0,1.2,2*nPML);

omegas = linspace(-pi*3e8/dx, pi*3e8/dx, 101);
zs = exp(-i*omegas*dt);

for depth = d
    
    k = kappaFn(depth);
    a = alphaFn(depth);
    sig = sigmaFn(depth);
    
    [numerator, denominator] = zCFSRIPML(dt, ...
        'kappa', k, ...
        'alpha', a,...
        'sigma', sig );
    
    sContinuous = sFn(k, sig, a, omegas);
    sFDTD = polyval(numerator, zs)./polyval(denominator, zs);
    
    figure(1)
    clf
    plot(omegas, real(sContinuous), '-', omegas, real(sFDTD), 'o');
    figure(2)
    clf
    plot(omegas, imag(sContinuous), '-', omegas, imag(sFDTD), 'o');
    pause
    
    
end

