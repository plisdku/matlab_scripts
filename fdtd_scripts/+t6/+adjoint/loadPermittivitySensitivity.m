function [dNumerExp, dDenomExp] = loadPermittivitySensitivity(freeDirection)
%function [dNumerExp, dDenomExp] = loadPermittivitySensitivity(freeDirection)

freeDirection = 1;

[numer1, denom1] = adjoint.readPermittivity('permittivity');
[dEpsilon, verts] = adjoint.readDeps();

% Extract the expected sensitivities into dNumerExp and dDenomExp

dNumerExp = 0*numer1;
dDenomExp = 0*denom1;

for xyz = 1:3
for vv = 1:length(verts)
    dEps = dEpsilon{verts(vv), freeDirection}.tensor{xyz,xyz};
    
    % dNumerExp comes from the EH part
    % dDenomExp comes from the DB part
    
    for lag = 1:length(dEps.EH)
        coeff = dEps.EH{lag}.coefficients;
        if ~isempty(coeff)
        pos = dEps.EH{lag}.positions + repmat([2 2 1], length(coeff),1);
        
        for pp = 1:length(coeff)
            dNumerExp(pos(pp,1), pos(pp,2), pos(pp,3), xyz, lag) = ...
                dNumerExp(pos(pp,1), pos(pp,2), pos(pp,3), xyz, lag) + ...
                coeff(pp);
        end
        end
    end
    
    for lag = 1:length(dEps.DB)
        coeff = dEps.DB{lag}.coefficients;
        
        if ~isempty(coeff)
        pos = dEps.DB{lag}.positions + repmat([2 2 1], length(coeff),1);
        
        for pp = 1:length(coeff)
            dDenomExp(pos(pp,1), pos(pp,2), pos(pp,3), xyz, lag) = ...
                dDenomExp(pos(pp,1), pos(pp,2), pos(pp,3), xyz, lag) + ...
                coeff(pp);
        end
        end
    end
end
end




