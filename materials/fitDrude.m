function [epsinf, omegap, gamma, errNorm] = fitDrude(epsilon, lambdas, varargin)
% [epsinf, omegap, gamma] = fitDrude(permittivity, lambdas) will find
% parameters for a Drude-model material that best approximate the tabulated
% data given.  These parameters will not be ideal for FDTD because they do
% not take into account the specific discretization of any simulation.
% [epsinf, omegap, gamma] = fitDrude(permittivity, lambdas, fitWeights)
% where fitWeights = [realWeight, imagWeight] lets you do manual "weight
% twiddling" on the fits, to improve the fit for real permittivity at the
% expense of the fit of imaginary permittivity, for instance.  Default
% value is [1 1], equal weight.

realWeight = 1;
imagWeight = 1;

if nargin > 2
    realWeight = varargin{1}(1);
    imagWeight = varargin{1}(2);
end

realImag = @(x) [realWeight*real(x); imagWeight*imag(x)];

% I have found that to getting a decent initial guess is essential.  The
% first guess is a systematic search over a wide range of possible epsinf,
% plasmon frequency and gamma.
epsinfs = linspace(-20, 20, 10);
wps = linspace(1e14, 1e16, 10);
gammas = linspace(1e14, 1e16, 10);

bestResnorm = inf;
for epsinf = epsinfs
for wp = wps
for gamma = gammas
    
    permGuess = drudePermittivity(lambdas, epsinf, wp, gamma);
    resnorm = norm(permGuess - epsilon);
    if resnorm < bestResnorm
        bestResnorm = resnorm;
        besteps = epsinf;
        bestwp = wp;
        bestgamma = gamma;
    end
    
end
end
end

% Second guess pass.  This uses the best guess from above to find the
% analytical fit (non-FDTD case).  Note that big numbers (wp and gamma) are
% scaled down to be around 1.  The epsinf, wp and gamma returned in this
% stage are what you would get for 

guess0 = [besteps bestwp bestgamma].*[1 1e-15 1e-15];

errEps1 = @(epsWpGamma, lambdas) ...
    epsilon - drudePermittivity(lambdas, epsWpGamma(1), epsWpGamma(2), ...
    epsWpGamma(3));

cplxErr1 = @(epsWpGamma) realImag(errEps1(epsWpGamma.*[1 1e15 1e15], lambdas));

opt = optimset('Display', 'off');

useMatlab = 1;
if useMatlab
    [guessParams, resnorm, residual] = lsqnonlin(cplxErr1, guess0, [], [], opt);
else
    lb = [-100 0 0];
    ub = [100 1e8 1e8];
    %lb = [epsinfs(1) wps(1) gammas(1)] .* [1 1e-15 1e-15];
    %ub = [epsinfs(end) wps(end) gammas(end)] .* [1 1e-15 1e-15];
    absErr = @(epsWpGamma) norm(cplxErr1(epsWpGamma));
    guessParams = pso.PSO(absErr, guess0, lb, ub);
end

errNorm = norm(cplxErr1(guessParams));

guessParams = guessParams.*[1, 1e15, 1e15];


epsinf = guessParams(1);
omegap = guessParams(2);
gamma = guessParams(3);
