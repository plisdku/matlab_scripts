function iii = bandwidth(x, varargin)
% bandwidth Return the frequency at which a spectrum drops below -3 dB.
%
%   bandwidth(x) returns the 3 dB point of the power abs(x).^2
%   bandwidth(x, dB) provides the level explicitly
%
% Example usage:
%   [srcSpectrum, freqs] = analysis.spectrum(srcFn(dt*(1:numTimesteps)), 'Dt', dt);
%   srcBandwidth = freqs(analysis.bandwidth(srcSpectrum(1:round(end/2))));
%   bandwidthLambda = 2*pi/srcBandwidth;
%
if nargin > 1
    dBlevel = varargin{1};
else
    dBlevel = -3;
end
powerLevel = 10^(dBlevel/10);

power = abs(x).^2;

maxPower = max(power);

iii = 2 + length(power) - ...
    find(power(end:-1:1) > maxPower*powerLevel, 1, 'first');