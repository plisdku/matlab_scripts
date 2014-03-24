function E_expression = multilayerSource(E_amplitudes, k_vecs, ...
    boundaries, var);

if nargin < 4
    var = 'z';
end

bnds = [-inf reshape(boundaries, 1, []) inf];

% In COMSOL we use exp(iwt - ikx).  Converting to the physics convention
% necessitates taking complex conjugates of everything.

layerExpr = @(ll) ...
    sprintf('(%s>%s[nm])*(%s<=%s[nm])*((%s)*exp(-i*((%s)[1/nm])*%s)+(%s)*exp(i*((%s)[1/nm])*%s))', ...
        var, num2str(bnds(ll)), var, num2str(bnds(ll+1)), ...
        num2str(conj(E_amplitudes(1,ll))), num2str(conj(k_vecs(ll))), var, ...
        num2str(conj(E_amplitudes(2,ll))), num2str(conj(k_vecs(ll))), var);

numLayers = numel(boundaries)+1;

E_expression = layerExpr(1);

for ll = 2:numLayers
    E_expression = [E_expression '+' layerExpr(ll)];
end




