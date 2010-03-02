function [numerator, denominator] = zCFSRIPML(dt, varargin)
% The permittivity of a material in FDTD may be written as a rational
% function in z space, where multiplication by z is a forward sample,
% division by z is a backward sample, and the time derivative by forward
% (backward) difference is written (z-1)/dt ( (1-1/z)/dt ).
%
% The product epsilon*mu, or permittivity*permeability, comes up again
% and again in FDTD analyses.  This function returns the product eps*mu
% for the Drude model, omitting the factor eps0*mu0 (this function will
% calculate RELATIVE permittivity and permeability).
%
% To evaluate this permittivity function, use polyval:
%
% permittivity = polyval(numerator, z)./polyval(denominator, z);

eps0 = 8.85e-12;
mu0 = 4e-7*pi;

% handle arguments
X.kappa = 1;
X.alpha = 0;
X.sigma = 0;

X = parseargs(X, varargin{:});
%X = parseargs(X, 'kappa', kappa, 'alpha', alpha, 'sigma', sigma);

numerator = ...
    [X.kappa*X.alpha + X.sigma + 2*X.kappa*eps0/dt, ...
    X.kappa*X.alpha + X.sigma - 2*X.kappa*eps0/dt ];

denominator = ...
    [X.alpha + 2*eps0/dt, ...
    X.alpha - 2*eps0/dt ];
