function [epsFDTD, epsAnalytical, freqs] = testPermittivity(varargin)

import t5.*;

X.Model = {'LossyDielectric', 'Dielectric', 'Drude'};
X.sigma = [];
X.epsr = [];
X.epsinf = [];
X.omegap = [];
X.tauc = [];
X.dxyz = [];
X.dt = [];
X.wavelengths = [600e-9 3000e-9];

X = parseargs(X, varargin{:});

validateattributes(X.dxyz, {'numeric'}, {'positive', 'vector'},...
    'testPermittivity', 'dxyz');
validateattributes(X.dt, {'numeric'}, {'positive', 'scalar'}, ...
    'testPermittivity', 'dt');

if numel(X.dxyz) ~= 3
    error('dxyz must be a three-element vector');
end

dxyz = X.dxyz;
dt = X.dt;
%dxyz = [5e-9 5e-9 5e-9];
%dt = 0.99*courant(dxyz);

wavelengthRange = [600e-9 3000e-9];
freqRange = 2*pi*3e8./wavelengthRange([2 1]);
centerFreq = mean(freqRange);
dFreq = diff(freqRange)/2;

t0 = 10/dFreq;
srcFn = @(t) exp(-(t-t0).^2*dFreq^2/4).*sin(centerFreq*(t-t0));

cutoffFreq = centerFreq + 5*dFreq; % a lot.

duration = 10*wavelengthRange(2)/3e8;
numT = round(duration/dt);

outputPeriod = floor(2*pi/cutoffFreq/dt);

gridRadius = duration*3e8/2;
gridRadiusYee = ceil(gridRadius/dxyz(1));

%%

%tt = linspace(0, duration, 10000);
%figure(1); clf
%plot(tt, srcFn(tt));

%[F, omega] = analysis.spectrum(srcFn(tt), 'Time', tt);

% figure(2); clf
% plot(omega, abs(F).^2);
% hold on
% plot([1 1]*cutoffFreq, [0, max(abs(F).^2)]);

% figure(3); clf
% plot(2*pi*3e8*1e9./omega, abs(F).^2);
% hold on
% plot([1 1]*2*pi*3e8*1e9/cutoffFreq, [0, max(abs(F).^2)])
% xlabel('wavelength (nm)')
% ylabel('F^2')
% xlim([0, wavelengthRange(2)*1e9])

%%

trogdor_begin(dxyz, dt, numT);

if strcmpi(X.Model, 'lossydielectric')
    validateattributes(X.epsr, {'numeric'}, {'scalar', 'positive', 'real'}, ...
        'testPermittivity', 'epsr');
    validateattributes(X.sigma, {'numeric'}, {'scalar', 'real'}, ...
        'testPermittivity', 'sigma');
    newLossyDielectric('material', 'sigma', X.sigma, 'epsr', X.epsr);
    
    permittivity = @(w) X.epsr + 1i*X.sigma./w/PhysicalConstants.eps0;
    
elseif strcmpi(X.Model, 'staticdielectric')
    validateattributes(X.epsr, {'numeric'}, {'scalar', 'positive', 'real'}, ...
        'testPermittivity', 'epsr');
    newDielectric('material', 'epsr', X.epsr);
    
    permittivity = @(w) X.epsr;
    
elseif strcmpi(X.Model, 'drude')
    validateattributes(X.epsinf, {'numeric'}, {'scalar', 'positive', 'real'}, ...
        'testPermittivity', 'epsinf');
    validateattributes(X.omegap, {'numeric'}, {'scalar', 'positive', 'real'}, ...
        'testPermittivity', 'omegap');
    validateattributes(X.tauc, {'numeric'}, {'scalar', 'positive', 'real'}, ...
        'testPermittivity', 'tauc');
    newDrude('material', 'epsinf', X.epsinf, 'omegap', X.omegap, ...
        'tauc', X.tauc);
    
    permittivity = @(w) X.epsinf - X.omegap^2./(w.^2 + 1i*w/X.tauc);
end

everywhere = [-gridRadiusYee, 0, 0, gridRadiusYee, 0, 0];

addGrid('Main');
addBlock('material', everywhere);

addCurrentSource('Field', 'jz', 'YeeCells', [0 0 0 0 0 0], 'TimeData', ...
    srcFn(dt*(0:numT-1)));

%addOutput('outAll', 'ez', 'YeeCells', everywhere);
addOutput('testData', 'ez', 'YeeCells', [1 0 0 2 0 0]);

trogdor_end('Directory', 'sim');

%%

!trogdor sim/params.xml > out.txt

%%

[F, freqs] = analysis.spectrum('testData');
F = squish(F);

delete('log.out.txt');
delete('out.txt');
delete('testData.txt');
delete('testData');

%%
wavevector = log(F(2,:)./F(1,:))/dxyz(1)/1i;
epsilon = wavevector.^2 ./ (PhysicalConstants.mu0 * freqs.^2);

epsFDTD = conj(epsilon / PhysicalConstants.eps0);
epsAnalytical = permittivity(freqs);




