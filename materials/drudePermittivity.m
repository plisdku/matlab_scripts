function permittivity = drudePermittivity(lambdas, epsinf, omegap, gamma)
% permittivity = drudePermittivity(lambdas, epsinf, omegap, gamma) will
% return the relative permittivity of a Drude-model material with the given
% asymptotic relative permittivity, plasmon frequency and relaxation
% constant gamma = 1/tauc.

omegas = 2*pi*2.99792458e8./lambdas;
tauc = 1/gamma;

permittivity = epsinf - omegap^2./(omegas.^2 + 1i*omegas/tauc);

