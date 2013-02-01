function [a w0] = bandLimitedGaussian(freqRange, power)
% bandLimitedGaussian  Return parameters for a modulated Gaussian function
% which falls to a given power at the edges of a given frequency range.
%
% [a w0] = bandLimitedGaussian(freqRange, power) returns the decay
% parameter a and center frequency w0 for a modulated Gaussian of the form
%
% sin(w0*t).*exp(-a*t.^2)
%

if power <= 0 || power >= 1
    error('Power must fall between 0 and 1');
end

dw = 0.5*(freqRange(2) - freqRange(1));

logP = log(power);

a = -0.5*dw^2/logP;
w0 = mean(freqRange);

% 
