function [numerator, denominator] = zDrude4(dt, varargin)
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


% handle arguments
X.omegap = 0;
X.tauc = 0;
X.epsinf = 0;

X = parseargs(X, varargin{:});
%X = parseargs(X, 'omegap', omegap, 'tauc', tauc, 'epsinf', epsinf);
X.omegap = X.omegap * sqrt(X.epsinf); % since Trogdor 4 is weird

% The numerator and denominator are both second-order polynomials.
epsilonNumerator = ...
    [X.epsinf*(2+dt/X.tauc),...
    2*X.omegap^2*dt^2 - 4*X.epsinf, ...
    X.epsinf*(2-dt/X.tauc)];
epsilonDenominator = ...
    [2+dt/X.tauc, ...
    -4, ...
    2-dt/X.tauc];

muNumerator = 1;
muDenominator = 1;

numerator = conv(epsilonNumerator, muNumerator);
denominator = conv(epsilonDenominator, muDenominator);



