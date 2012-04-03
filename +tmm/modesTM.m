function [k_better, t22, kvals] = modesTM(varargin)
% modesTM  Find the wavevectors of bound TM modes of a multilayer
% dielectric structure, using the reflection pole method (RPM).
%
% k = modesTM(boundaries, epsr, mur, omega, kMin, kMax) returns an array
% of complex wavevectors, each corresponding to a likely bound mode of the
% system.
%
%   k               array of interface-parallel wavevectors [1/m]
%   boundaries      z positions of boundaries between materials [m]
%   epsr            relative permittivity of each layer [unitless]
%   mur             relative permeability of each layer [unitless]
%   omega           optical angular frequency [1/s]
%   kMin, kMax      bounds on real parts of k to search over.  A good
%                   choice for kMin is k0/nMax, where k0 is the free-space
%                   wavevector of light (omega/c) and nMax is the larger of
%                   the indices of refraction of the first and last layers
%                   of the multilayer. [1/m]
%
% [k, t22, kvals] = modesTM(...) returns the internal data used to locate
% the waveguide modes.  The reflection pole method works by searching for
% zeros of the (2,2) component of a transfer matrix by inspecting the phase
% of this element (t22) while varying the parallel wavevector from kMin to
% kMax.
%
%   t22             the transfer matrix element T(2,2) at each point in the
%                   wavevector sweep [unitless]
%   kvals           the wavevectors searched by the method
%
%
% The candidate mode wavevectors may be used with solveTM() to obtain field
% profiles for each mode.  This example extracts the field profile for the
% first detected mode.
%
% EXAMPLE:
%
% omega = 2*pi*3e8/600e-9;
% k0 = omega/3e8;
% boundaries = [0e-9];
% epsr = [1, -10 + 1i]; % a dielectric-metal interface
% mur = [1 1];
% kMin = 1.01*k0;
% kMax = 4*k0;
% outPositions = linspace(-1000e-9, 1000e-9);
%
% kParallels = tmm.modesTM(boundaries, epsr, mur, omega, kMin, kMax);
% [hx ey ez] = tmm.solveTM(boundaries, epsr, mur, omega, kParallels(1), ...
%       outPositions, true);
%
% plot(outPositions*1e9, real(hx), outPositions*1e9, imag(hx));
% xlabel('z (nm)');
% ylabel('H_x (nm)');
% legend('Real', 'Imag')
%

warning('modesTM is deprecated, please use modes()')

[k_better, t22, kvals] = tmm.modes(varargin{:}, 'tm');