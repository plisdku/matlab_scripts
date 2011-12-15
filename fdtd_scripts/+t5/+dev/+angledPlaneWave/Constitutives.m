classdef Constitutives
    
    methods (Static)
        
        function func = drude(epsinf, omegap, tauc, dt)
            func = @(z) PhysicalConstants.eps0 * (epsinf + ...
                2*dt^2*omegap^2./((2+dt/tauc)-(2-dt/tauc)./z)./(z-1));
        end
        
        function func = dielectric(epsr)
            func = @(z) PhysicalConstants.eps0*epsr*ones(size(z));
        end
        
        function func = conductor(epsr, sigma, dt)
            func = @(z) PhysicalConstants.eps0*epsr + dt*sigma./(z-1);
        end
        
        function func = permeability(mur)
            func = @(z) PhysicalConstants.mu0*mur*ones(size(z));
        end
    end
    
end