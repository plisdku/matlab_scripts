function y = evalMaterial(numer, denom, z)

y = polyval(numer, z)./polyval(denom, z);
