function [numerator, denominator] = materialModel(modelName)

dielectricNumer = @(epsRel) epsRel;
dielectricDenom = @(epsRel) 1;

conductorNumer = @(epsInf, loss) [epsInf+loss, loss-epsInf];
conductorDenom = @(epsInf, loss) [1, -1]; % maybe it'll work?

debyeNumer = @(epsS, epsInf, tau, dt) ...
    [(epsS*dt+2*tau*epsInf)/(dt+2*tau), (epsS*dt-2*tau*epsInf)/(dt+2*tau)];
debyeDenom = @(epsS, epsInf, tau, dt) [1, (dt-2*tau)/(dt+2*tau)];

drudeNumer = @(epsInf, omegap, gamma, dt) ...
    [epsInf*(2+gamma*dt), -4*epsInf + 2*omegap^2*dt^2, epsInf*(2-gamma*dt)];
drudeDenom = @(epsInf, omegap, gamma, dt) ...
    [2+gamma*dt, -4, 2-gamma*dt];

lorentzNumer = @(epsinf, epss, omegap, gamma, dt) ...
    [ omegap^2*dt^2*epss + 2*gamma*dt*epsinf + 2, ...
      -4*epsinf, ...
      omegap^2*dt^2*epss - 2*gamma*dt*epsinf + 2];
lorentzDenom = @(epsinf, epss, omegap, gamma, dt) ...
    [ omegap^2*dt^2 + 2*gamma*dt + 2, ...
      -4, ...
      omegap^2*dt^2 - 2*gamma*dt + 2];
  
switch lower(modelName)
    case 'debye'
        numerator = debyeNumer;
        denominator = debyeDenom;
    case 'conductor'
        numerator = conductorNumer;
        denominator = conductorDenom;
    case 'drude'
        numerator = drudeNumer;
        denominator = drudeDenom;
    case 'lorentz'
        numerator = lorentzNumer;
        denominator = lorentzDenom;
    case 'dielectric'
        numerator = dielectricNumer;
        denominator = dielectricDenom;
    otherwise
        error(sprintf(['Material model must be debye, conductor,', ...
            'drude, dielectric or lorentz (received %s'], modelName));
end
