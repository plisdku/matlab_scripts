function [srcE, srcH] = generate(direction, omega, dielectricEpsilon, modeIndex,...
    srcFn)
% generate     Create TM waveguide mode source for multilayer stack.
% Usage: generate(direction, frequency, dielectricEpsilon)
%        generate(direction, frequency, dielectricEpsilon, modeIndex)
%        generate(direction, frequency, dielectricEpsilon, modeIndex, srcFn)
%

import t5.dev.multilayerMode.*

tfsfreq_0;

posE = {afp.posEx, afp.posEy, afp.posEz};
posH = {afp.posHx, afp.posHy, afp.posHz};

[afpE afpH] = phasors(direction, omega, dielectricEpsilon, 1, posE, posH);

lambda = 2*pi*3e8./omega;
k0 = 2*pi/lambda;
kMin = 1.01*k0;
kMax = 3*k0;

% I think providing a default srcFn is a bad idea, now.
%if ~exist('srcFn', 'var')
%    t0 = 8*lambda/3e8/afp.dt;
%    width = 4*lambda/3e8/afp.dt;
%    
%    srcFn = @(t) 0.5*(1+erf((t-t0)/width)).*exp(-1i*t*omega);
%end

%% Transform to the time domain and write to file

fprintf('Writing fields to file.\n');
srcE = @(xyz,n) real(srcFn(afp.dt*n) .* afpE{xyz});
srcH = @(xyz,n) real(srcFn(afp.dt*(n+0.5)) .* afpH{xyz});

fh = fopen(afp.inputFile, 'w');

for tt = 1:afp.numT
    
    for xyz = 1:3
        fwrite(fh, srcH(xyz,tt), 'float32');
    end
    
    for xyz = 1:3
        fwrite(fh, srcE(xyz,tt), 'float32');
    end
    
end
fclose(fh);

