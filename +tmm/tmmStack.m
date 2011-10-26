function [E, H, T, R, epsrOut, murOut] = tmmStack(boundaries, epsr, mur,...
    inputE, omega, ks, outputPosE, outputPosH)
% Usage:
% [E, H, T, R, epsrOut] = tmmStack(boundaries, epsr, mur, inputE, omega, ks,
%   outputPosE, outputPosH)
%
% E is an array of E fields measured at outputPosE
% H is an array of H fields measured at outputPosH
% T is the transmitted power from 0 to 1
% R is the reflected power from 0 to 1
% epsrOut is an array of epsr values measured at outputPosE
% 
% boundaries is an array of positions where E and H are continuous [meters]
% 
% ks is an array of wavenumbers, one per layer, including the media before
% and after the multilayer. [1/meters]
% 
% eps is an array of relative permittivities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% mu is an array of relative permeabilities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% inputE is either the amplitude of incoming E at the left, or a
% two-element array with the amplitudes of incoming E at the left and
% right, respectively.  PRESENTLY ONLY FORWARD AMPLITUDE IS SUPPORTED.
% [arbitrary units]
% 
% outputPosE is a vector of positions to evaluate the E field at [meters]
% 
% outputPosH is a vector of positions to evaluate the H field at [meters]

import tmm.*;

%boundaries = [0e-9, 100e-9, 150e-9];
%epsr = [1, 5, 10, 1];
%mur = [1, 1, 1, 1];
n = sqrt(epsr.*mur);
%inputE = 1;
%omega = 2*pi*3e8/30e-9;
%ks = n*omega/3e8;

%outputPosE = linspace(-50e-9, 200e-9, 500);
%outputPosH = outputPosE;

% Make matrices to convert [E+, E-] in layer n to [E+, E] in layer n+1,
% and vice-versa

% E(n+1) = forwardMatrix{n}*E(n)
forwardMatrix = cell(length(boundaries), 1);

for nLayer = 1:length(boundaries)
    forwardMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundaries(nLayer)) * ...
        matrixEE2EH(omega, ks(nLayer), boundaries(nLayer));
end

%% Get incoming and reflected fields consistent with unit transmission
% amplitude

% E_last = transferEE * E_first
transferEE = eye(2);

% E_N = transferLayer{N} * E_first
transferLayer = cell(length(boundaries)+1, 1); % for intermediate layers

transferLayer{1} = transferEE;
for nLayer = 1:length(forwardMatrix)
    transferEE = forwardMatrix{nLayer}*transferEE;
    transferLayer{nLayer+1} = transferEE;
end

inverseTransferEE = inv(transferEE); % E_first = inverseTransferEE * E_last

E0 = inverseTransferEE*[1;0];
t = 1/E0(1);
r = E0(2) / E0(1);
T = abs(t)^2;
R = abs(r)^2;
assert( abs( 1.0 - T - R ) < 1e-6 ); % t^2 + r^2 = 1

%% Get the forward and backward E in each layer (use transferLayer)

E = 0*outputPosE;
H = 0*outputPosH;
epsrOut = E;
murOut = H;

intervals = [-inf, boundaries(1,:), inf];
for nLayer = 1:length(boundaries)+1
    indicesE = find( outputPosE > intervals(nLayer) & ...
        outputPosE <= intervals(nLayer+1));
    indicesH = find(outputPosH > intervals(nLayer) & ...
        outputPosH <= intervals(nLayer+1));
    
    En = transferLayer{nLayer}*E0/E0(1)*inputE;
    
    for ii = indicesE
        z = outputPosE(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z);
        EH = ee2eh*En;
        E(ii) = EH(1);
    end
    
    for ii = indicesH
        z = outputPosH(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z);
        EH = ee2eh*En;
        H(ii) = EH(2);
    end
    
    epsrOut(indicesE) = epsr(nLayer);
    murOut(indicesH) = mur(nLayer);
end




