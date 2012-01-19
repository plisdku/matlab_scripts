% Test the multilayer solver

import tmm.*;

mu0 = 4e-7*pi;
eps0 = 8.854187817e-12;
c = 1/sqrt(eps0*mu0);

checkSmall = @(a) assert(abs(a) < 1e-5);
checkRelativelyClose = @(a, b) checkSmall(norm(a-b)/norm(a));
checkClose = @(a, b) checkSmall(norm(a-b));

lambda = 1000e-9;
k = 2*pi/lambda;
omega = c*k;


%% Test reflection from an interface

boundaries = [0];
epsr = [1,3];
mur = [1, 1];
n = sqrt(epsr.*mur);
inputE = 1;

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
    
    [Ex, Hy, Hz, T, R] = solveTE(boundaries, epsr, mur, inputE, omega, ...
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
    
    [Hx, Ey, Ez, T, R] = solveTM(boundaries, epsr, mur, omega, ...
        kParallel);
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
inputE = 1;

lambdas = linspace(450e-9, 900e-9, 100);
reflections = 0*lambdas;

for ii = 1:length(lambdas)
    omega = 2*pi*c/lambdas(ii);
    [Ex, Hy, Hz, T, R] = solveTE(boundaries, epsr, mur, inputE, omega, 0);
    reflections(ii) = R;
end

reflection650 = interp1(lambdas, reflections, 650e-9);
checkSmall(reflection650);

disp('Quarter-wave antireflection coating: passed');


























