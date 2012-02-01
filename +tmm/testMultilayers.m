% Test the multilayer solver

import tmm.*;

mu0 = 4e-7*pi;
eps0 = 8.854187817e-12;
eta0 = sqrt(mu0/eps0);
c = 1/sqrt(eps0*mu0);

checkSmall = @(a) assert(abs(a) < 1e-5);
checkRelativelyClose = @(a, b) checkSmall(norm(a-b)/norm(a));
checkRelativelyClose1 = @(a,b) checkSmall(norm(a-b,1)/norm(a,1)); % use 1-norm
checkClose = @(a, b) checkSmall(norm(a-b));

lambda = 1000e-9;
k = 2*pi/lambda;
omega = c*k;


%% Test reflection from an interface

boundaries = [0];
epsr = [1,3];
mur = [1, 1];
n = sqrt(epsr.*mur);

ks = 2*pi*n/lambda;

fresnelR_s = @(n1, n2, theta1, theta2) ...
    (n1*cos(theta1) - n2*cos(theta2))^2 /...
    (n1*cos(theta1) + n2*cos(theta2))^2;
fresnelR_p = @(n1, n2, theta1, theta2) ...
    (n1*cos(theta2) - n2*cos(theta1))^2 /...
    (n1*cos(theta2) + n2*cos(theta1))^2;

thetas = linspace(0,pi/2 - 0.001,100);

zs = linspace(-4e-6, 4e-6, 1000);

%% TE

reflections = 0*thetas;
for nn = 1:length(thetas)
    theta1 = thetas(nn);
    theta2 = asin(n(1)*sin(theta1)/n(2));
    checkClose(n(1)*sin(theta1), n(2)*sin(theta2)); % check Snell's law
    
    kParallel = sin(theta1)*ks(1);
    
    [Ex, Hy, Hz, T, R] = solveTE(boundaries, epsr, mur, omega, ...
        kParallel);
    reflections(nn) = R;
    checkRelativelyClose(R, fresnelR_s(n(1), n(2), theta1, theta2));
end

disp('Fresnel equations, s-polarization: passed');

%% TM

reflections = 0*thetas;
for nn = 1:length(thetas)
    theta1 = thetas(nn);
    theta2 = asin(n(1)*sin(theta1)/n(2));
    checkClose(n(1)*sin(theta1), n(2)*sin(theta2)); % check Snell's law
    
    kParallel = sin(theta1)*ks(1);
    
    [Hx, Ey, Ez, T, R] = solveTM(boundaries, epsr, mur, omega, kParallel);
    reflections(nn) = R;
    checkRelativelyClose(R, fresnelR_p(n(1), n(2), theta1, theta2));
end

disp('Fresnel equations, p-polarization: passed');

%% Test antireflection layer

n1 = 1;
n3 = 3;
n2 = sqrt(n1*n3);
n = [n1 n2 n3];
epsr = n.^2;
mur = [1 1 1];

lambda = 650e-9;

boundaries = [0 lambda/4/n2];

lambdas = linspace(450e-9, 900e-9, 100);
reflections = 0*lambdas;

for ii = 1:length(lambdas)
    omega = 2*pi*c/lambdas(ii);
    [Ex, Hy, Hz, T, R] = solveTE(boundaries, epsr, mur, omega, 0);
    reflections(ii) = R;
end

reflection650 = interp1(lambdas, reflections, 650e-9);
checkSmall(reflection650);

disp('Quarter-wave antireflection coating: passed');

%% Test Maxwell's equations for TE and TM fields
% I launch a wave with a complex k-vector against a multilayer with a lossy
% material in it and check that Maxwell's equations hold, i.e. that the
% E fields are compatible with the H field for TM incidence, and that the H
% fields are compatible with the E field for TE incidence.

n1 = 1;
n2 = 3 + 1i;
n3 = 2;
n = [n1 n2 n3];
epsr = n.^2;
mur = [1 1 1];
lambda = 500e-9;
omega = 2*pi*c/lambda;
k0 = 2*pi/lambda;
ky = k0/3 * exp(1i*0.01);
boundaries = [0 100e-9];
outBounds = [-100e-9, boundaries, 200e-9];

compPlot = @(y1, y2) plot(zPos, y1, '-', zPos, y2, '-');
err2Plot = @(y1, y2) plot(zPos, abs(y1-y2).^2);
errHist = @(y1, y2) hist(abs(y1-y2));

norm1 = @(y1, y2) norm(y1-y2, 1)/norm(y1,1);
norm2 = @(y1, y2) norm(y1-y2, 2)/norm(y1,2);

% I do a pair of checks of Maxwell's equations here.
% One of them uses a numerical derivative of a transverse field component
% along the layer-normal direction.  The gradient spikes at boundaries
% between materials and is generally less accurate at the endpoints of its
% domain.  I get around this issue by checking Maxwell's equations in each
% layer separately, first of all, and secondarily by using a 1-norm to
% evaluate the error instead of a 2-norm.  The basic message is the same
% (the error is small) but the 1-norm makes the measurement more robust to
% outliers.

for layer = 1:numel(n)
    zPos = linspace(outBounds(layer)+1e-9, outBounds(layer+1)-1e-9, 1000);
    [Ex, Hy, Hz, T, R, epsX, muY, muZ] = solveTE(...
        boundaries, epsr, mur, omega, ky, zPos);
    
    %fprintf('TE Maxwell check, layer %i: ', layer);
    checkRelativelyClose(-omega*mu0*muZ.*Hz, ky*Ex);
    checkRelativelyClose1(-1i*omega*mu0*muY.*Hy, -gradient(Ex, zPos));
    %fprintf('passed\n');
end

disp('TE Maxwell equations: passed');

for layer = 1:numel(n)
    zPos = linspace(outBounds(layer)+1e-9, outBounds(layer+1)-1e-9, 1000);
    [Hx, Ey, Ez, T, R, muX, epsY, epsZ] = solveTM(...
        boundaries, epsr, mur, omega, ky, zPos);
    
    %fprintf('TM Maxwell check, layer %i: ', layer);
    checkRelativelyClose(-omega*eps0*epsZ.*Ez, -ky*Hx);
    checkRelativelyClose1(-1i*omega*eps0*epsY.*Ey, gradient(Hx, zPos));
    %fprintf('passed\n');
end

disp('TM Maxwell equations: passed');








