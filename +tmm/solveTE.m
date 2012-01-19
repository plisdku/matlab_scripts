function [Ex, Hy, Hz, T, R, epsrEx, murHy, murHz, transferEE] = solveTE(boundaries, ...
    epsr, mur, inputE, omega, kParallel, varargin)
% Usage:
% [Ex, Hy, Hz, T, R, epsrEx, murHy, murHz, transferMatrix] = solveTE(boundaries, epsr, mur,
%   inputE, omega, ky, outputPosEx, outputPosHy, outputPosHz)
%
% Ex is an array of transverse E fields measured at outputPosEx
% Hy is an array of transverse H fields measured at outputPosHy
% Hz is an array of normal H fields measured at outputPosHz
% T is the transmitted power from 0 to 1
% R is the reflected power from 0 to 1
% epsrEx is an array of epsr values measured at outputPosEx
% murHy is an array of mur values measured at outputPosHy
% murHz is an array of mur values measured at outputPosHz
% 
% boundaries is an array of positions where E and H are continuous [meters]
% 
% ky is the k vector parallel to the boundary. [1/meters]
%
% epsr is an array of relative permittivities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% mur is an array of relative permeabilities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% inputE is the amplitude of incoming E at the left [arbitrary units]
% 
% outputPosEx is a vector of positions to evaluate the Ex field at.
% (optional) [meters]
% 
% outputPosHy is a vector of positions to evaluate the Hy field at.
% (optional) [meters]
%
% outputPosHz is a vector of positions to evaluate the Hz field at.
% (optional) [meters]

import tmm.*;

outputPosEx = [];
outputPosHy = [];
outputPosHz = [];

if nargin > 6
    outputPosEx = reshape(varargin{1},1,[]);
    if nargin > 7
        outputPosHy = reshape(varargin{2},1,[]);
        if nargin > 8
            outputPosHz = reshape(varargin{3},1,[]);
        end
    end
end

mu0 = 4e-7*pi;
eps0 = 8.854187817e-12;
c = 1/sqrt(eps0*mu0);

n = sqrt(epsr.*mur);

ks = sqrt(omega^2*n.^2/c^2 - kParallel^2);
ks(imag(ks) < 0) = -ks(imag(ks) < 0); % decay goes the right way now

% Make matrices to convert [E+, E-] in layer n to [E+, E-] in layer n+1,
% and vice-versa

% E(n+1) = forwardMatrix{n}*E(n)
forwardMatrix = cell(length(boundaries), 1);

for nLayer = 1:length(boundaries)
    forwardMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundaries(nLayer), mur(nLayer+1)) * ...
        matrixEE2EH(omega, ks(nLayer), boundaries(nLayer), mur(nLayer));
end

%% Get incoming and reflected fields consistent with unit transmission

% E_last = transferEE * E_first
transferEE = eye(2);

% E_N = transferLayer{N} * E_first
transferLayer = cell(length(boundaries)+1, 1); % for intermediate layers

transferLayer{1} = transferEE;
for nLayer = 1:length(forwardMatrix)
    transferEE = forwardMatrix{nLayer}*transferEE;
    transferLayer{nLayer+1} = transferEE;
end

% The transfer matrix is transferEE.

E0 = transferEE \ [1;0];
t = 1/E0(1);
r = E0(2) / E0(1);
T = abs(t)^2;
R = abs(r)^2;
%assert( abs( 1.0 - T - R ) < 1e-6 ); % t^2 + r^2 = 1 % not generally true

%% Get the forward and backward E in each layer (use transferLayer)

Ex = 0*outputPosEx;
Hy = 0*outputPosHy;
Hz = 0*outputPosHz;
epsrEx = Ex;
murHy = Hy;
murHz = Hz;

intervals = [-inf, boundaries(:)', inf];
for nLayer = 1:length(boundaries)+1
    
    indicesEx = [];
    indicesHy = [];
    indicesHz = [];
    
    if ~isempty(outputPosEx)
        indicesEx = find( outputPosEx > intervals(nLayer) & ...
            outputPosEx <= intervals(nLayer+1));
    end
    
    if ~isempty(outputPosHy)
        indicesHy = find(outputPosHy > intervals(nLayer) & ...
         outputPosHy <= intervals(nLayer+1));
    end
    
    if ~isempty(outputPosHz)
        indicesHz = find(outputPosHz > intervals(nLayer) & ...
         outputPosHz <= intervals(nLayer+1));
    end
    
    En = transferLayer{nLayer}*E0/E0(1)*inputE;
    
    for ii = indicesEx
        z = outputPosEx(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Ex(ii) = EH(1);
    end
    
    for ii = indicesHy
        z = outputPosHy(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Hy(ii) = EH(2);
    end
    
    for ii = indicesHz
        z = outputPosHz(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Hz(ii) = -EH(1)*kParallel/omega/mur(nLayer)/mu0;
        %Hz(ii) = EH(2)*kParallel/ks(nLayer); % WRONG
    end
    
    epsrEx(indicesEx) = epsr(nLayer);
    murHy(indicesHy) = mur(nLayer);
    murHz(indicesHz) = mur(nLayer);
end


