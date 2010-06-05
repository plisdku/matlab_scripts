classdef PhysicalConstants
    
    properties (Constant = true)
        eps0 = 8.854187817e-12;
        mu0 = 4e-7*pi;
        c = 1/sqrt(PhysicalConstants.eps0*PhysicalConstants.mu0);
        eta0 = sqrt(PhysicalConstants.mu0 / PhysicalConstants.eps0);
        hbar = 1.05457148e-34;
        boltzmann = 1.3806503e-23;
        eV = 1.60217646e-19; % Joules
        electronMass = 9.10938188e-31;
        avogadro = 6.0221415e23;
    end
    
end