function addPMLParams(grid, gridXML, doc)

%{
    allDirectionsDefault["sigma"] =
        "(d^3)*0.8*4/(((mu0/eps0)^0.5)*dx)";
    allDirectionsDefault["alpha"] =
        "d*3e8*eps0/(50*dx)";
    allDirectionsDefault["kappa"] =
        "1 + (5-1)*(d^3)";
%}

params = grid.PMLParams;

if ~isempty(params.alpha) || ~isempty(params.sigma) || ...
    ~isempty(params.kappa)
    
    paramsXML = doc.createElement('PML');
    
    if ~isempty(params.alpha)
        paramsXML.setAttribute('alpha', params.alpha);
    end
    
    if ~isempty(params.kappa)
        paramsXML.setAttribute('kappa', params.kappa);
    end
    
    if ~isempty(params.sigma)
        paramsXML.setAttribute('sigma', params.sigma);
    end
    
    gridXML.appendChild(paramsXML);
end

