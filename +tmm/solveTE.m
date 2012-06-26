function [Ex, Hy, Hz, T, R, epsrEx, murHy, murHz, transferEE] = solveTE(...
    boundary_z, epsr, mur, omega, kParallel, varargin)
% Usage:
% [Ex, Hy, Hz, T, R, epsrEx, murHy, murHz, transferMatrix] = solveTE(boundary_z, epsr, mur,
%   omega, ky, output_z, forceBoundModes)
%
% Ex is an array of transverse E fields measured at output_zEx
% Hy is an array of transverse H fields measured at output_zHy
% Hz is an array of normal H fields measured at output_zHz
% T is the transmitted power from 0 to 1
% R is the reflected power from 0 to 1
% epsrEx is an array of epsr values measured at output_zEx
% murHy is an array of mur values measured at output_zHy
% murHz is an array of mur values measured at output_zHz
% transferMatrix is the 2x2 transfer matrix mapping left- and
%   right-propagating E amplitudes on the left side of the stack to left- and
%   right-propagating E amplitudes on the right side of the stack
% 
% boundary_z is an array of z positions where E and H are continuous [meters]
% 
% ky is the k vector parallel to the boundary. [1/meters]
%
% epsr is an array of relative permittivities, one per layer, including the
% media before and after the multilayer.  Positive imaginary permittivity 
% connotes loss.  [unitless]
% 
% mur is an array of relative permeabilities, one per layer, including the
% media before and after the multilayer [unitless]
% 
% output_z is a vector of z positions to evaluate the fields at.  It may
% also be a cell array with three vectors of positions, one each for Ex, Hy
% and Hz. (optional) [meters]
%
% forceBoundModes can be true or false (false by default).  If true, the
% transfer matrices and field amplitudes will be adjusted so no inbound
% waves are present. (optional)
%
% A forward-propagating wave is represented as exp(1i*(k*k - w*t)).  This
% negative frequency convention implies that lossy materials must have
% positive imaginary permittivities.

import tmm.*;

output_z = [];
output_zEx = [];
output_zHy = [];
output_zHz = [];

if numel(varargin) > 0
    output_z = varargin{1};
    if iscell(output_z)
        output_zEx = output_z{1};
        output_zHy = output_z{2};
        output_zHz = output_z{3};
    else
        output_zEx = output_z;
        output_zHy = output_z;
        output_zHz = output_z;
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

% Make matrices to convert [E+, E-] in layer n to [E+, E-] in layer n+1,
% and vice-versa

% E(n+1) = forwardMatrix{n}*E(n)
forwardMatrix = cell(length(boundary_z), 1);

for nLayer = 1:length(boundary_z)
    forwardMatrix{nLayer} = ...
        matrixEH2EE(omega, ks(nLayer+1), boundary_z(nLayer), mur(nLayer+1)) * ...
        matrixEE2EH(omega, ks(nLayer), boundary_z(nLayer), mur(nLayer));
end

%% Get incoming and reflected fields consistent with unit transmission

% E_last = transferEE * E_first
transferEE = eye(2);

% E_N = transferLayer{N} * E_first
transferLayer = cell(length(boundary_z)+1, 1); % for intermediate layers

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

%% Rescale to unit incidence

E0 = E0 / E0(1);

%% For mode solutions:

normalizationPos = [];

if forceBoundModes
    E0(1) = 0;
    transferLayer{end}(2,:) = 0;
    
    normalizationPos = linspace(boundary_z(1) - 3*2*pi/abs(ks(1)), ...
        boundary_z(end) + 3*2*pi/abs(ks(end)), 10000);
end

%% Get the forward and backward E in each layer (use transferLayer)

Ex = 0*output_zEx;
Hy = 0*output_zHy;
Hz = 0*output_zHz;
epsrEx = Ex;
murHy = Hy;
murHz = Hz;

layerEnergy = [];

intervals = [-inf, boundary_z(:)', inf];
for nLayer = 1:length(boundary_z)+1
    
    indicesEx = [];
    indicesHy = [];
    indicesHz = [];
    indicesNormalization = [];
    
    if ~isempty(output_zEx)
        indicesEx = find( output_zEx > intervals(nLayer) & ...
            output_zEx <= intervals(nLayer+1));
    end
    
    if ~isempty(output_zHy)
        indicesHy = find(output_zHy > intervals(nLayer) & ...
         output_zHy <= intervals(nLayer+1));
    end
    
    if ~isempty(output_zHz)
        indicesHz = find(output_zHz > intervals(nLayer) & ...
         output_zHz <= intervals(nLayer+1));
    end
    
    if ~isempty(normalizationPos)
        indicesNormalization = find(normalizationPos > intervals(nLayer) & ...
            normalizationPos <= intervals(nLayer+1));
    end
    
    En = transferLayer{nLayer}*E0;
    
    for ii = indicesEx
        z = output_zEx(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Ex(ii) = EH(1);
    end
    
    for ii = indicesHy
        z = output_zHy(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Hy(ii) = EH(2);
    end
    
    for ii = indicesHz
        z = output_zHz(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Hz(ii) = -EH(1)*kParallel/omega/mur(nLayer)/mu0;
    end
    
    for ii = indicesNormalization
        z = normalizationPos(ii);
        ee2eh = matrixEE2EH(omega, ks(nLayer), z, mur(nLayer));
        EH = ee2eh*En;
        Ex_normalization(ii) = EH(1);
        Hz_normalization(ii) = -EH(1)*kParallel/omega/mur(nLayer)/mu0;
    end
    
    epsrEx(indicesEx) = epsr(nLayer);
    murHy(indicesHy) = mur(nLayer);
    murHz(indicesHz) = mur(nLayer);
    
    %{
    if forceBoundModes
        %Hn = matrixEE2HH_TE(omega, ks(nLayer), mur(nLayer)) * En;
        Hn = -En * kParallel/omega/mur(nLayer)/mu0;
        
        z = normalizationPos(indicesNormalization);
        figure(20); clf
        hLayer = Hn(1)*exp(1i*ks(nLayer)*z) + Hn(2)*exp(-1i*ks(nLayer)*z);
        title('H')
        plot(z, real(hLayer), z, imag(hLayer))
        figure(21); clf
        eLayer = En(1)*exp(1i*ks(nLayer)*z) + En(2)*exp(-1i*ks(nLayer)*z);
        plot(z, real(eLayer), z, imag(eLayer));
        title('E')
        
        fprintf('En    %2.4g    %2.4g\n', En(1), En(2));
        fprintf('Hn    %2.4g    %2.4g\n', Hn(1), Hn(2));
        
        l1 = intervals(nLayer);
        l2 = intervals(nLayer+1);
        
        if nLayer == 1 % first layer: only include outward waves
            layerEnergy(nLayer) = 1/(2i*ks(nLayer)) * En(2)*Hn(2) * ( ...
                exp(-2i*ks(nLayer)*l2));
        elseif nLayer == length(boundary_z)+1
            layerEnergy(nLayer) = 1/(2i*ks(nLayer)) * En(1)*Hn(1) * ( ...
                -exp(2i*ks(nLayer)*l1));
        else
            term(1) = (l2-l1)*(En(1)*Hn(2) + En(2)*Hn(1));
            term(2) = 1/(2i*ks(nLayer)) * En(1)*Hn(1)*( ...
                exp(2i*ks(nLayer)*l2) - exp(2i*ks(nLayer)*l1));
            term(3) = 1/(2i*ks(nLayer)) * En(2)*Hn(2)*( ...
                exp(-2i*ks(nLayer)*l2) - exp(-2i*ks(nLayer)*l1));
            layerEnergy(nLayer) = sum(term);
        end
    end
    %}
end

%%


if ~isempty(normalizationPos)
    modeEnergy = trapz(normalizationPos, -Hz_normalization.*Ex_normalization);

    Ex = Ex / abs(sqrt(modeEnergy));
    Hy = Hy / abs(sqrt(modeEnergy));
    Hz = Hz / abs(sqrt(modeEnergy));
end

%fprintf('Mode energy %2.4g     Integrated %2.4g\n', modeEnergy, sum(layerEnergy));
%keyboard
