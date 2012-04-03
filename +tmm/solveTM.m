function [Hx, Ey, Ez, T, R, murHx, epsrEy, epsrEz, transferHH] = solveTM(...
    boundary_z, epsr, mur, omega, kParallel, varargin)
% Usage:
% [Hx, Ey, Ez, T, R, murHx, epsrEy, epsrEz, transferMatrix] = solveTM(boundary_z, epsr,
%   mur, omega, ky, output_z, forceBoundModes, normalizationPos)
%
% Hx is an array of transverse H fields measured at output_zHx
% Ey is an array of transverse E fields measured at output_zEy
% Ez is an array of normal E fields measured at output_zEz
% T is the transmitted power from 0 to 1
% R is the reflected power from 0 to 1
% murHx is an array of mur values measured at output_zHx
% epsrEy is an array of epsr values measured at output_zEy
% epsrEz is an array of epsr values measured at output_zEz
% transferMatrix is the 2x2 transfer matrix mapping left- and
%   right-propagating H amplitudes on the left side of the stack to left- and
%   right-propagating H amplitudes on the right side of the stack
% 
% boundary_z is an array of z positions where E and H are continuous [meters]
% 
% ky is the k vector parallel to the boundary. [1/meters]
% 
% epsr is an array of relative permittivities, one per layer, including the
% media before and after the multilayer.  Positive imaginary permittivity 
% connotes loss. [unitless]
% 
% mur is an array of relative permeabilities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% output_z is a vector of z positions to evaluate the fields at.  It may
% also be a cell array with three vectors of positions, one each for Hx, Ey
% and Ez. (optional) [meters]
%
% forceBoundModes can be true or false (false by default).  If true, the
% transfer matrices and field amplitudes will be adjusted so no inbound
% waves are present.  The modes returned will be normalized to unit "power".
% (optional)
%
% A forward-propagating wave is represented as exp(1i*(k*k - w*t)).  This
% negative frequency convention implies that lossy materials must have
% positive imaginary permittivities.

import tmm.*;

output_z = [];
output_zHx = [];
output_zEy = [];
output_zEz = [];

if numel(varargin) > 0
    output_z = varargin{1};
    if iscell(output_z)
        output_zHx = reshape(output_z{1}, 1, []);
        output_zEy = reshape(output_z{2}, 1, []);
        output_zEz = reshape(output_z{3}, 1, []);
    else
        output_zHx = reshape(output_z, 1, []);
        output_zEy = reshape(output_z, 1, []);
        output_zEz = reshape(output_z, 1, []);
    end
end

forceBoundModes = false;
if numel(varargin) > 1
    forceBoundModes = varargin{2};
end

mu0 = 4e-7*pi;
eps0 = 8.854187817e-12;
c = 1/sqrt(eps0*mu0);

n = sqrt(epsr.*mur);

ks = sqrt(omega^2*n.^2/c^2 - kParallel^2);
ks(imag(ks) < 0) = -ks(imag(ks) < 0); % decay goes the right way now

% Make matrices to convert [H+, H-] in layer n to [H+, H-] in layer n+1,
% and vice-versa

% H(n+1) = forwardMatrix{n}*H(n)
forwardMatrix = cell(length(boundary_z), 1);

for nLayer = 1:length(boundary_z)
    forwardMatrix{nLayer} = ...
        matrixHE2HH(omega, ks(nLayer+1), boundary_z(nLayer), epsr(nLayer+1))* ...
        matrixHH2HE(omega, ks(nLayer), boundary_z(nLayer), epsr(nLayer));
end

%% Get incoming and reflected fields consistent with unit transmission

% H_last = transferHH * H_first
transferHH = eye(2);

% H_N = transferLayer{N} * H_first
transferLayer = cell(length(boundary_z)+1, 1); % for intermediate layers

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

%% Rescale to unit incidence

H0 = H0 / H0(1);

%% For mode solutions:

normalizationPos = [];

if forceBoundModes
    H0(1) = 0;
    transferLayer{end}(2,:) = 0;
    
    normalizationPos = linspace(boundary_z(1) - 3*2*pi/abs(ks(1)), ...
        boundary_z(end) + 3*2*pi/abs(ks(end)), 10000);
end

%% Get the forward and backward E in each layer (use transferLayer)

Hx = 0*output_zHx;
Ey = 0*output_zEy;
Ez = 0*output_zEz;
murHx = Hx;
epsrEy = Ey;
epsrEz = Ez;

Ez_normalization = 0*normalizationPos;
Hx_normalization = 0*normalizationPos;

intervals = [-inf, boundary_z(:)', inf];
for nLayer = 1:length(boundary_z)+1
    
    indicesHx = [];
    indicesEy = [];
    indicesEz = [];
    indicesNormalization = [];
    
    if ~isempty(output_zHx)
        indicesHx = find( output_zHx > intervals(nLayer) & ...
            output_zHx <= intervals(nLayer+1));
    end
    
    if ~isempty(output_zEy)
        indicesEy = find(output_zEy > intervals(nLayer) & ...
         output_zEy <= intervals(nLayer+1));
    end
    
    if ~isempty(output_zEz)
        indicesEz = find(output_zEz > intervals(nLayer) & ...
         output_zEz <= intervals(nLayer+1));
    end
    
    if ~isempty(normalizationPos)
        indicesNormalization = find(normalizationPos > intervals(nLayer) & ...
            normalizationPos <= intervals(nLayer+1));
    end
    
    Hn = transferLayer{nLayer}*H0;
    
    for ii = indicesHx
        z = output_zHx(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Hx(ii) = HE(1);
    end
    
    for ii = indicesEy
        z = output_zEy(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Ey(ii) = HE(2);
    end
    
    for ii = indicesEz
        z = output_zEz(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Ez(ii) = HE(1)*kParallel/omega/epsr(nLayer)/eps0;
    end
    
    for ii = indicesNormalization
        z = normalizationPos(ii);
        hh2he = matrixHH2HE(omega, ks(nLayer), z, epsr(nLayer));
        HE = hh2he*Hn;
        Hx_normalization(ii) = HE(1);
        Ez_normalization(ii) = HE(1)*kParallel/omega/epsr(nLayer)/eps0;
    end
    
    murHx(indicesHx) = mur(nLayer);
    epsrEy(indicesEy) = epsr(nLayer);
    epsrEz(indicesEz) = epsr(nLayer);
end


if ~isempty(normalizationPos)
    modeEnergy = trapz(normalizationPos, Ez_normalization.*Hx_normalization);

    Hx = Hx / sqrt(modeEnergy);
    Ey = Ey / sqrt(modeEnergy);
    Ez = Ez / sqrt(modeEnergy);
end



