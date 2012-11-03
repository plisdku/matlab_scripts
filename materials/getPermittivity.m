function [eps, lam] = getPermittivity(matName, varargin)
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

lamRange = [-inf, inf];
if nargin >= 2
    lamRange = varargin{1};
    assert(length(lamRange) == 2);
end

[lam, ~, ~, er, ei] = loadData(matName);

[lam, indices] = sort(lam);
er = er(indices);
ei = ei(indices);
eps = er + 1i*ei;

indices = find(lam >= lamRange(1) & lam <= lamRange(2));

lam = lam(indices);
eps = eps(indices);


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