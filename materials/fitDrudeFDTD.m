function [epsinf, omegap, gamma] = fitDrudeFDTD(epsilon, lambdas, dxdydzdt, varargin)
% [epsinf, omegap, gamma] = fitDrudeFDTD(permittivity, lambdas, dxdydzdt) will
% attempt to find material parameters for a simulated Drude model that best
% approximates the given complex permittivity.  The parameters will not be
% the same as the best fit in the "real world" due to discretization.
% [epsinf, omegap, gamma] = fitDrudeFDTD(permittivity, lambdas, dxdydzdt,
% ndims) will try to adjust the parameters to minimize grid dispersion as
% well in ndims dimensions.  The default value is one.
% [epsinf, omegap, gamma] = fitDrudeFDTD(permittivity, lambdas, dxdydzdt,
% ndims, fitWeights) where fitWeights = [realWeight, imagWeight] lets you do 
% manual "weight twiddling" on the fits, to improve the fit for real
% permittivity at the expense of the fit of imaginary permittivity, for 
% instance.  Default value is [1 1], equal weight.

ndims = 1;
realWeight = 1;
imagWeight = 1;

if nargin > 3
    ndims = varargin{1};
end

if nargin > 4
    realWeight = varargin{2}(1);
    imagWeight = varargin{2}(2);
end

%epsilon = eps;
%lambdas = lam;
%dxdydzdt = [dx dx dx dt];

% We need a guess to start with.  Get one from fitDrude.  It's pretty good.

[epsinf, omegap, gamma] = fitDrude(epsilon, lambdas, [realWeight, imagWeight])

% Now do the final fit: correct for the discretization.  In my experience,
% the medium-scale least-squares algorithm (option 'LargeScale', 'off')
% converges a lot faster than the large-scale algorithm.

realImag = @(x) [realWeight*real(x); imagWeight*imag(x)];

if ndims == 1
    errEps = @(epsWpGamma, lambdas, dxdydzdt) ...
        epsilon ...
        - drudePermittivityFDTD(lambdas, [1 0 0], dxdydzdt, epsWpGamma(1), ...
        epsWpGamma(2), epsWpGamma(3));
elseif ndims == 2
    errEps = @(epsWpGamma, lambdas, dxdydzdt) ...
        [epsilon; epsilon] - ...
        [   drudePermittivityFDTD(lambdas, [1 0 0], dxdydzdt, epsWpGamma(1), ...
            epsWpGamma(2), epsWpGamma(3)); ...
            drudePermittivityFDTD(lambdas, [1 1 0], dxdydzdt, epsWpGamma(1), ...
            epsWpGamma(2), epsWpGamma(3)) ];
elseif ndims == 3
    errEps = @(epsWpGamma, lambdas, dxdydzdt) ...
        [epsilon; epsilon; epsilon] - ...
        [   drudePermittivityFDTD(lambdas, [1 0 0], dxdydzdt, epsWpGamma(1), ...
            epsWpGamma(2), epsWpGamma(3)); ...
            drudePermittivityFDTD(lambdas, [1 1 0], dxdydzdt, epsWpGamma(1), ...
            epsWpGamma(2), epsWpGamma(3)); ...
            drudePermittivityFDTD(lambdas, [1 1 1], dxdydzdt, epsWpGamma(1), ...
            epsWpGamma(2), epsWpGamma(3)) ];
end
    

cplxErr = @(epsWpGamma) realImag(errEps(epsWpGamma.*[1,1e15,1e15],...
    lambdas, dxdydzdt));

initialGuess = [epsinf, omegap, gamma].*[1, 1e-15, 1e-15];
opt = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt');
[theParams, resnorm, residual] = lsqnonlin(cplxErr, initialGuess, [], [], opt);
%%

epsinf = theParams(1);
omegap = theParams(2)*1e15;
gamma = theParams(3)*1e15;

