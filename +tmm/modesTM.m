function [k_better, t22, kvals] = modesTM(boundaries, epsr, mur, omega, kMin, kMax)

k0 = omega/3e8;

%% Sweep the k-vector and collect reflection denominator terms

kvals = linspace(kMin, kMax, 10000);
nEff = kvals/k0;

t22 = 0*kvals;
for nn = 1:length(kvals)
    t22(nn) = reflectionDenominator(boundaries, epsr, mur, omega, kvals(nn));
end

%%

phiRPM = unwrap(angle(t22));

%% Make a list of plausible resonances

yy = gradient(phiRPM, nEff);
yy = yy / max(yy(:));

inds = findLocalMax(yy, 3);
if inds(1) == 1
    inds = inds(2:end);
end
if inds(end) == length(yy)
    inds = inds(1:end-1);
end

candidateHeights = yy(inds);
candidateCenters = nEff(inds);
candidateWidths = 0*candidateCenters;

for nn = 1:length(inds)
    [iLeft, iRight] = findFWHM(yy, inds(nn));
    candidateWidths(nn) = 0.5*(spline(1:length(nEff), nEff, iRight) - ...
        spline(1:length(nEff), nEff, iLeft));
end

%figure(1010)
%plot(nEff, yy);
%keyboard

%% Fit with Lorentzians
% lorentzians() is a subroutine

amplitudes = @(p) p(1:length(p)/3);
centers = @(p) p(length(p)/3+1:2*length(p)/3);
gammas = @(p) p(2*length(p)/3+1:end);
testFn = @(p,x)lorentzians(x, amplitudes(p), centers(p), gammas(p));

numBumps = length(candidateCenters);

pMin = [0*candidateHeights, candidateCenters-0.1, 0*candidateWidths];
pMax = [3*candidateHeights, candidateCenters+0.1, 1];

pInit = [candidateHeights, candidateCenters, candidateWidths];

%myPlotFn = @(x, optimValues, state) durrplot(1010, nEff, yy, x, pInit, testFn);
%opts = optimset('PlotFcns', myPlotFn);

soln = struct;
[soln.pFinal, soln.resnorm, soln.residual, soln.exitflag, soln.output] =...
    lsqcurvefit(testFn, pInit, nEff, yy, pMin, pMax);

fits.amplitudes = amplitudes(soln.pFinal);
fits.centers = centers(soln.pFinal);
fits.gammas = gammas(soln.pFinal);


%figure(1)
%clf
%plot(nEff, yy, '-', ...
%    nEff, testFn(soln.pFinal, nEff), 'o', ...
%    nEff, testFn(pInit, nEff), '--');
%legend('t22', 'Final', 'Guess')


%% Improve the fits
%
% This step is actually pretty important.  Sometimes the Lorentzian fits do
% not provide sufficiently accurate estimates of the zeros of t22.  I can
% zoom in locally and use a minimization algorithm to improve the accuracy.
% Using this technique I've had some success in forcing asymmetrical mode
% profiles to look symmetrical when they're supposed to be symmetrical.

n_fits = fits.centers + 1i*fits.gammas;
k_fits = n_fits*k0;

t22_func = @(k) reflectionDenominator(boundaries, epsr, mur, omega, k);

for nn = 1:length(n_fits)
    
    kScale = [real(k_fits(nn)), imag(k_fits(nn))];
    kStart = [real(k_fits(nn)), imag(k_fits(nn))];
    
    kScaledStart = kStart ./ kScale;
    kScaledStart(isnan(kScaledStart)) = 0;
    
    t22_vec = @(k2) reflectionDenominator(boundaries, epsr, mur, omega, ...
        k2(1)*kScale(1) + 1i*k2(2)*kScale(2));
    
    %minMe = @(k) abs(t22_func(k));
    minMe = @(kScaled) abs(t22_vec(kScaled));
    %opt = optimset('TolFun', 1e-11, 'TolX', k0*1e-8);
    %n_betters(nn) = fminsearch(minMe, k_fits(nn))/k0;
    k2Better = fminsearch(minMe, kScaledStart);
    
    k_better(nn) = k2Better(1)*kScale(1) + 1i*k2Better(2)*kScale(2);
    n_better(nn) = k_better(nn)/k0;
end



function y = lorentzians(x, amplitudes, x0, gamma)

% this is the usual definition
%lor = @(x, x0, gamma) (1/pi)*gamma./( (x-x0).^2 + gamma^2);

lor = @(x, x0, gamma) gamma^2./( (x-x0).^2 + gamma^2);

nL = length(amplitudes);

y = amplitudes(1)*lor(x, x0(1), gamma(1));

for nn = 2:nL
    y = y + amplitudes(nn)*lor(x, x0(nn), gamma(nn));
end


function t22 = reflectionDenominator(boundaries, epsr, mur, omega, k)

[hx, ey, ez, T, R, murHx, epsrEy, epsrEz, txMatrix] = ...
    tmm.solveTM(boundaries, epsr, mur, omega, k);

t22 = txMatrix(2,2);

    %tmm.solveTM('Boundaries', boundaries, ...
    %    'Permittivity', epsr, 'Permeability', mur, ...
    %    'Frequency', omega, 'Wavevector', k);
    %tmm.solveTM(boundaries, epsr, mur, 1, omega, k);

