% Test the transfer matrix stuff

import tmm.*;

checkSmall = @(a) assert(abs(a) < 1e-5);
checkRelativelyClose = @(a, b) checkSmall(norm(a-b)/norm(a));
checkClose = @(a, b) checkSmall(norm(a-b));

lambda = 1000;
k = 2*pi/lambda;
omega = k;
x = 100;

%% Test basic matrices

% TE, mur = 1

ee2eh = matrixEE2EH(omega, k, x);
eh2ee = matrixEH2EE(omega, k, x);
checkRelativelyClose(ee2eh, inv(eh2ee));

propag = matrixPropagateEH(omega, k, x);
propagTheoretical = matrixEE2EH(omega, k, x)*matrixEH2EE(omega, k, 0);
checkRelativelyClose(propag, propagTheoretical);

% TE, mur = 3.21 (arbitrary)

mur = 3.21;

ee2eh = matrixEE2EH(omega, k, x, mur);
eh2ee = matrixEH2EE(omega, k, x, mur);
checkRelativelyClose(ee2eh, inv(eh2ee));

propag = matrixPropagateEH(omega, k, x, mur);
propagTheoretical = matrixEE2EH(omega, k, x, mur) * ...
    matrixEH2EE(omega, k, 0, mur);
checkRelativelyClose(propag, propagTheoretical);

disp('EE to EH conversion: passed');

% TM, epsr = 1

hh2he = matrixHH2HE(omega, k, x);
he2hh = matrixHE2HH(omega, k, x);
checkRelativelyClose(hh2he, inv(he2hh));

propag = matrixPropagateHE(omega, k, x);
propagTheoretical = matrixHH2HE(omega, k, x)*matrixHE2HH(omega, k, 0);
checkRelativelyClose(propag, propagTheoretical);

% TM, epsr = 4.19 (arbitrary)

epsr = 4.19;

hh2he = matrixHH2HE(omega, k, x, epsr);
he2hh = matrixHE2HH(omega, k, x, epsr);
checkRelativelyClose(hh2he, inv(he2hh));

propag = matrixPropagateHE(omega, k, x, epsr);
propagTheoretical = matrixHH2HE(omega, k, x, epsr) * ...
    matrixHE2HH(omega, k, 0, epsr);
checkRelativelyClose(propag, propagTheoretical);

disp('HH to HE conversion: passed');

%% Test that transfer matrices are unimodular

checkRelativelyClose(det(matrixPropagateEH(omega,k,x)), 1);
checkRelativelyClose(det(matrixPropagateHE(omega,k,x)), 1);
checkRelativelyClose(det(matrixPropagateEH(omega,k,x, 1.4)), 1);
checkRelativelyClose(det(matrixPropagateHE(omega,k,x, 2.2)), 1);

disp('Unimodular transfer matrices: passed');

%% Test that propagation is invertible

checkRelativelyClose(matrixPropagateEH(omega,k,-x),...
    inv(matrixPropagateEH(omega,k,x)));
checkRelativelyClose(matrixPropagateHE(omega,k,-x),...
    inv(matrixPropagateHE(omega,k,x)));

checkRelativelyClose(matrixPropagateEH(omega,k,-x,1.1),...
    inv(matrixPropagateEH(omega,k,x,1.1)));
checkRelativelyClose(matrixPropagateHE(omega,k,-x,44),...
    inv(matrixPropagateHE(omega,k,x,44)));

disp('Forward-backward propagation: passed');

%% Test phase accumulation

mur = 4;
E = [1; 0]; % forward wave
EH = matrixPropagateEH(omega,k,x,mur)*matrixEE2EH(omega,k,0,mur)*E;
checkRelativelyClose(EH(1), exp(1i*k*x));

epsr = 3;
H = [0; 1]; % backward wave
HE = matrixPropagateHE(omega,k,x,epsr)*matrixHH2HE(omega,k,0,epsr)*H;
checkRelativelyClose(HE(1), exp(-1i*k*x));

disp('Phase due to propagation: passed');

%% Check that TE and TM give same results at normal incidence

epsr = 2 + .3i;
mur = 1.1 + 2i;
n = sqrt(epsr*mur);
eta = sqrt(mur/epsr); % impedance of medium

k = omega*n; % k vector in the medium, at normal incidence

propagateTM = matrixPropagateHE(omega, k, x, epsr);
propagateTE = matrixPropagateEH(omega, k, x, mur);

% Matrix to swap E to H and H to -E, the similarity relation for TE and TM.
swapEH = [0 -1; 1 0];

% I'll test both directions because I'm superstitious.  (-:
checkRelativelyClose(propagateTM, swapEH*propagateTE*swapEH');
checkRelativelyClose(swapEH*propagateTM*swapEH', propagateTE);

disp('TE and TM agreement at normal incidence: passed');

%% Check that TE and TM give the Fresnel equations

fresnelR_s = @(n1, n2, theta1, theta2) ...
    (n1*cos(theta1) - n2*cos(theta2))^2 /...
    (n1*cos(theta1) + n2*cos(theta2))^2;
fresnelR_p = @(n1, n2, theta1, theta2) ...
    (n1*cos(theta2) - n2*cos(theta1))^2 /...
    (n1*cos(theta2) + n2*cos(theta1))^2;

epsr1 = 1;
epsr2 = 5;
mur = 1;
n1 = sqrt(epsr1);
n2 = sqrt(epsr2);

lambda = 1000;
k = 2*pi/lambda;
omega = k;

thetas = linspace(0, 0.49*pi, 100);

refTE = 0*thetas;
refTM = 0*thetas;
for ii = 1:length(thetas)
    theta1 = thetas(ii);
    theta2 = asin(n1*sin(theta1)/n2);
    
    k1 = omega*n1;
    k2 = omega*n2;
    
    checkClose(n1*sin(theta1), n2*sin(theta2));
    
    fresnelMatrixTE = ...
        matrixEH2EE(omega, k1*cos(theta1), 0, mur) * ...
        matrixEE2EH(omega, k2*cos(theta2), 0, mur);
    fresnelMatrixTM = ...
        matrixHE2HH(omega, k1*cos(theta1), 0, epsr1) * ...
        matrixHH2HE(omega, k2*cos(theta2), 0, epsr2);
    
    % [E+, E-] (layer 1) = fresnelMatrixTE * [E+, E-] (layer 2)
    % [H+, H-] (layer 1) = fresnelMatrixTM * [H+, H-] (layer 2)
    
    refTE(ii) = fresnelMatrixTE(2,1) / fresnelMatrixTE(1,1);
    refTM(ii) = fresnelMatrixTM(2,1) / fresnelMatrixTM(1,1);
    
    checkRelativelyClose(refTE(ii)^2, fresnelR_s(n1, n2, theta1, theta2));
    checkRelativelyClose(refTM(ii)^2, fresnelR_p(n1, n2, theta1, theta2));
    
end

disp('TE and TM Fresnel equations: passed');

