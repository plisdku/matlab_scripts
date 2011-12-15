function iii = bandwidth(x, varargin)
% bandwidth Return the frequency at which a spectrum drops below -3 dB
%
%   bandwidth(x) returns the 3 dB point
%   bandwidth(x, dB) provides the level explicitly

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