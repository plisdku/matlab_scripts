function [Hx, Ey, Ez, T, R, murHx, epsrEy, epsrEz, transferHH] = solveTM(...
    boundaries, epsr, mur, omega, kParallel, varargin)
% Usage:
% [Hx, Ey, Ez, T, R, murHx, epsrEy, epsrEz, transferMatrix] = solveTM(boundaries, epsr,
%   mur, omega, kParallel, outputPos, forceBoundModes)
%
% Hx is an array of transverse H fields measured at outputPosHx
% Ey is an array of transverse E fields measured at outputPosEy
% Ez is an array of normal E fields measured at outputPosEz
% T is the transmitted power from 0 to 1
% R is the reflected power from 0 to 1
% murHx is an array of mur values measured at outputPosHx
% epsrEy is an array of epsr values measured at outputPosEy
% epsrEz is an array of epsr values measured at outputPosEz
% 
% boundaries is an array of positions where E and H are continuous [meters]
% 
% kParallel is the k vector parallel to the boundary. [1/meters]
% 
% epsr is an array of relative permittivities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% mur is an array of relative permeabilities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% outputPos is a vector of positions to evaluate the fields at.  It may
% also be a cell array with three vectors of positions, one each for Hx, Ey
% and Ez. (optional) [meters]
%
% forceBoundModes can be true or false (false by default).  If true, the
% transfer matrices and field amplitudes will be adjusted so no inbound
% waves are present. (optional)

outputPos = [];
outputPosHx = [];
outputPosEy = [];
outputPosEz = [];

if numel(varargin) > 0
    outputPos = varargin{1};
    if iscell(outputPos)
        outputPosHx = outputPos{1};
        outputPosEy = outputPos{2};
        outputPosEz = outputPos{3};
    else
        outputPosHx = outputPos;
        outputPosEy = outputPos;
        outputPosEz = outputPos;
    end
end

forceBoundModes = false;
if numel(varargin) > 1
    forceBoundModes = varargin{2};
end

import tmm.*;

mu0 = 4e-7*pi;
eps0 = 8.854187817e-12;
c = 1/sqrt(eps0*mu0);

n = sqrt(epsr.*mur);

ks = sqrt(omega^2*n.^2/c^2 - kParallel^2);
ks(imag(ks) < 0) = -ks(imag(ks) < 0); % decay goes the right way now

% Make matrices to convert [H+, H-] in layer n to [H+, H-] in layer n+1,
% and vice-versa

% H(n+1) = forwardMatrix{n}*H(n)
forwardMatrix = cell(length(boundaries), 1);

for nLayer = 1:length(boundaries)
    forwardMatrix{nLayer} = ...
        matrixHE2HH(omega, ks(nLayer+1), boundaries(nLayer), epsr(nLayer+1))* ...
        matrixHH2HE(omega, ks(nLayer), boundaries(nLayer), epsr(nLayer));
end

%% Get incoming and reflected fields consistent with unit transmission

% H_last = transferHH * H_first
transferHH = eye(2);

% H_N = transferLayer{N} * H_first
transferLayer = cell(length(boundaries)+1, 1); % for intermediate layers

transferLayer{1} = transferHH;
for nLayer = 1:length(forwardMatrix)
    transferHH = forwardMatrix{nLayer}*transferHH;
    transferLayer{nLayer+1} = transferHH;
end

H0 = transferHH \ [1;0];
t = 1/H0(1);
r = H0(2) / H0(1);
T = abs(t)^2;
R = abs(r)^2;
%assert( abs( 1.0 - T - R ) < 1e-6 ); % t^2 + r^2 = 1 % not generally true

%% For mode solutions:

%if X.ForceBoundMode
if forceBoundModes
    H0(1) = 0;
    transferLayer{end}(2,:) = 0;
end

%% Get the forward and backward E in each layer (use transferLayer)

Hx = 0*outputPosHx;
Ey = 0*outputPosEy;
Ez = 0*outputPosEz;
murHx = Hx;
epsrEy = Ey;
epsrEz = Ez;

intervals = [-inf, boundaries(:)', inf];
for nLayer = 1:length(boundaries)+1
    
    indicesHx = [];
    indicesEy = [];
    indicesEz = [];
    
    if ~isempty(outputPosHx)
        indicesHx = find( outputPosHx > intervals(nLayer) & ...
            outputPosHx <= intervals(nLayer+1));
    end
    
    if ~isempty(outputPosEy)
        indicesEy = find(outputPosEy > intervals(nLayer) & ...
         outputPosEy <= intervals(nLayer+1));
    end
    
    if ~isempty(outputPosEz)
        indicesEz = find(outputPosEz > intervals(nLayer) & ...
         outputPosEz <= intervals(nLayer+1));
    end
    
    Hn = transferLayer{nLayer}*H0;
    
    for ii = indicesHx
        z = outputPosHx(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Hx(ii) = HE(1);
    end
    
    for ii = indicesEy
        z = outputPosEy(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Ey(ii) = HE(2);
    end
    
    for ii = indicesEz
        z = outputPosEz(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Ez(ii) = -HE(2)*kParallel/ks(nLayer);
    end
    
    murHx(indicesHx) = mur(nLayer);
    epsrHy(indicesEy) = epsr(nLayer);
    epsrHz(indicesEz) = epsr(nLayer);
end




