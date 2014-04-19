function [eps, lambdas] = getPermittivity(matName, varargin)
% [epsilon, lambda] = getPermittivity(materialName) will return tabulated
% real and imaginary relative permittivity over the full available range of
% wavelengths.  Epsilon is in units of eps0 = 8.85e-12 and lambda is in
% meters.
% [epsilon, lambda] = getPermittivity(materialName, lambdaRange) where
% lambdaRange = [lambdaMin, lambdaMax] will truncate the returned data to a
% given wavelength range.
%
% Available materials include
%   Ag
%   Al
%   Al32Cr
%   Au
%   AuGa
%   CNT
%   Co_hex_vis+IR
%   Cr
%   diamond
%   GST-A
%   GST-C
%   SiO2
%   Ti
%   TiO2

    X.WavelengthRange = [-inf, inf];
    X.Wavelengths = [];

    if nargin == 2
        X.WavelengthRange = varargin{1};
    else
        X = parseargs(X, varargin{:});
    end
    
    [lambdas, ~, ~, er, ei] = loadData(matName);

    [lambdas, indices] = sort(lambdas);
    eps = er(indices) + 1i*ei(indices);

    if ~isempty(X.Wavelengths)
        eps = spline(lambdas, eps, X.Wavelengths);
    else
        if length(X.WavelengthRange) ~= 2
            error('Wavelength range must be a 2-element vector.');
        end
        
        ii = find(lambdas >= X.WavelengthRange(1) & ...
            lambdas <= X.WavelengthRange(2));
        lambdas = lambdas(ii);
        eps = eps(ii);
    end

end


function [lam n k er ei] = loadData(matName)

    thisfile = mfilename('fullpath');
    [pathstr, name] = fileparts(thisfile);

    fname = [pathstr, filesep, 'txt', filesep, matName, '.txt'];

    try
        AA = dlmread(fname);
    catch
        error('Cannot find data file for material named %s', matName);
    end

    lam = AA(:,1)*1e-9;
    n = AA(:,2);
    k = AA(:,3);
    er = AA(:,4);
    ei = AA(:,5);

end