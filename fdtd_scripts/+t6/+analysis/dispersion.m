function func = dispersion(numer, denom)
% dispersion  Return a function f(w) = polyval(numer, 1i*w)./polyval(denom, 1i*w)
%
% Example: static conductor with eps = 3.3 - 5.5i at 800 nm wavelength
%   numer = [3.3 0.0432];
%   denom = [1 0];
%   epsilon = dispersion(numer, denom);
%
%   epsilon(2*pi/800)
%
%       ans =
%
%          3.3000 - 5.5000i
% 
% Note that under the positive frequency convention a lossy material has a
% negative imaginary permittivity.

func = @(w) polyval(numer, 1i*w) ./ polyval(denom, 1i*w);